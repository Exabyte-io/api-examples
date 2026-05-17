"""
Mat3ra Prompt-to-Code Server
=============================
A tiny HTTP server that accepts a natural-language prompt and returns
the equivalent Python code that calls the Mat3ra API.

POST /ask   { "prompt": "List my silicon materials" }
            -> { "code": "...", "tool": "list_materials", "args": {...} }

GET  /      -> simple HTML playground

Run:
    python mcp_server/code_server.py

Then try:
    curl -s -X POST http://localhost:8888/ask \\
         -H "Content-Type: application/json" \\
         -d '{"prompt": "List my silicon materials"}' | python3 -m json.tool
"""

from __future__ import annotations

import json
import re
import textwrap
from http.server import BaseHTTPRequestHandler, HTTPServer

# ---------------------------------------------------------------------------
# Intent detection + code generation (rule-based, no LLM needed for the PoC)
# ---------------------------------------------------------------------------

BOILERPLATE = """\
import os
from mat3ra.api_client.endpoints.{module} import {cls}

ACCOUNT_ID = os.environ["MAT3RA_ACCOUNT_ID"]
AUTH_TOKEN  = os.environ["MAT3RA_AUTH_TOKEN"]
OWNER_ID    = os.environ.get("MAT3RA_ORGANIZATION_ID") or ACCOUNT_ID

endpoint = {cls}(
    "platform.mat3ra.com", 443,
    ACCOUNT_ID, AUTH_TOKEN,
    "2018-10-01", secure=True
)
"""


def detect_intent(prompt: str) -> dict:
    """
    Map a natural-language prompt to a { tool, args, module, cls, call } dict.
    In a production system this would be done by an LLM; here we use simple
    keyword patterns so the server works with zero external dependencies.
    """
    p = prompt.lower()

    # ── materials ────────────────────────────────────────────────────────────
    # list / search materials
    formula_match = re.search(
        r"\b([A-Z][a-z]?(?:\d*[A-Z][a-z]?)*\d*)\b", prompt
    )
    formula = formula_match.group(1) if formula_match else None

    if re.search(r"\b(list|search|find|get|show)\b.*\bmaterial", p):
        args: dict = {"owner._id": "OWNER_ID"}
        if formula and formula not in ("List", "Show", "Find", "Get"):
            args["formula"] = formula
        query_str = json.dumps(args, indent=4).replace('"OWNER_ID"', "OWNER_ID")
        call = f"materials = endpoint.list(\n    {query_str}\n)"
        return dict(
            tool="list_materials", args=args,
            module="materials", cls="MaterialEndpoints", call=call,
        )

    if re.search(r"\b(create|add|upload|make)\b.*\bmaterial", p):
        call = textwrap.dedent("""\
            config = {
                "name": "My Material",
                "basis": {
                    "elements": [{"id": 1, "value": "Si"}, {"id": 2, "value": "Si"}],
                    "coordinates": [{"id": 1, "value": [0, 0, 0]},
                                    {"id": 2, "value": [0.25, 0.25, 0.25]}],
                    "units": "crystal",
                },
                "lattice": {"type": "FCC", "a": 3.867, "b": 3.867, "c": 3.867,
                            "alpha": 60, "beta": 60, "gamma": 60,
                            "units": {"length": "angstrom", "angle": "degree"}},
                "tags": ["MCP"],
            }
            material = endpoint.create(config, OWNER_ID)""")
        return dict(
            tool="create_material", args={"config": "..."},
            module="materials", cls="MaterialEndpoints", call=call,
        )

    if re.search(r"\bget\b.*\bmaterial\b.*\bid\b|\bmaterial\b.*\bid\b.*\bget\b", p):
        id_match = re.search(r"\b([a-f0-9]{24})\b", prompt)
        mid = id_match.group(1) if id_match else "<material_id>"
        call = f'material = endpoint.get("{mid}")'
        return dict(
            tool="get_material", args={"material_id": mid},
            module="materials", cls="MaterialEndpoints", call=call,
        )

    # ── workflows ────────────────────────────────────────────────────────────
    if re.search(r"\b(list|search|find|show)\b.*\bworkflow", p):
        call = 'workflows = endpoint.list({"owner._id": OWNER_ID})'
        return dict(
            tool="list_workflows", args={},
            module="workflows", cls="WorkflowEndpoints", call=call,
        )

    # ── jobs ─────────────────────────────────────────────────────────────────
    if re.search(r"\b(list|search|find|show)\b.*\bjob", p):
        call = 'jobs = endpoint.list({"owner._id": OWNER_ID})'
        return dict(
            tool="list_jobs", args={},
            module="jobs", cls="JobEndpoints", call=call,
        )

    if re.search(r"\b(create|submit|run|launch)\b.*\bjob", p):
        call = textwrap.dedent("""\
            config = {
                "owner":     {"_id": OWNER_ID},
                "_material": {"_id": "<material_id>"},
                "workflow":  {"_id": "<workflow_id>"},
                "name":      "MCP Job",
            }
            job = endpoint.create(config)
            endpoint.submit(job["_id"])
            job = endpoint.get(job["_id"])   # re-fetch to see status""")
        return dict(
            tool="create_and_submit_job",
            args={"material_id": "<material_id>", "workflow_id": "<workflow_id>"},
            module="jobs", cls="JobEndpoints", call=call,
        )

    if re.search(r"\b(get|status|check)\b.*\bjob", p):
        id_match = re.search(r"\b([a-f0-9]{24})\b", prompt)
        jid = id_match.group(1) if id_match else "<job_id>"
        call = f'job = endpoint.get("{jid}")'
        return dict(
            tool="get_job", args={"job_id": jid},
            module="jobs", cls="JobEndpoints", call=call,
        )

    # ── fallback ─────────────────────────────────────────────────────────────
    return dict(
        tool="unknown", args={},
        module="materials", cls="MaterialEndpoints",
        call="# Could not map prompt to a Mat3ra API call.\n"
             "# Supported: list/create/get materials, list workflows, list/create/get jobs.",
    )


def build_code(intent: dict) -> str:
    header = BOILERPLATE.format(module=intent["module"], cls=intent["cls"])
    return header.rstrip() + "\n\n" + intent["call"] + "\n"


# ---------------------------------------------------------------------------
# HTTP server
# ---------------------------------------------------------------------------

HTML_PLAYGROUND = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Mat3ra Prompt-to-Code</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { background: #0d1117; color: #c9d1d9; font-family: system-ui, sans-serif;
         display: flex; flex-direction: column; align-items: center;
         min-height: 100vh; padding: 40px 16px; }
  h1 { font-size: 1.6rem; margin-bottom: 8px; color: #58a6ff; }
  p  { color: #8b949e; margin-bottom: 24px; font-size: 0.9rem; }
  .card { background: #161b22; border: 1px solid #30363d; border-radius: 12px;
          padding: 24px; width: 100%; max-width: 720px; }
  textarea { width: 100%; background: #0d1117; color: #c9d1d9;
             border: 1px solid #30363d; border-radius: 8px;
             padding: 12px; font-size: 0.95rem; resize: vertical;
             outline: none; min-height: 72px; font-family: inherit; }
  textarea:focus { border-color: #58a6ff; }
  button { margin-top: 12px; background: #238636; color: #fff;
           border: none; border-radius: 8px; padding: 10px 24px;
           font-size: 0.95rem; cursor: pointer; transition: background 0.2s; }
  button:hover { background: #2ea043; }
  .label { font-size: 0.75rem; color: #8b949e; margin: 16px 0 6px; }
  pre  { background: #0d1117; border: 1px solid #30363d; border-radius: 8px;
         padding: 16px; font-size: 0.85rem; overflow-x: auto;
         white-space: pre-wrap; color: #79c0ff; min-height: 48px; }
  .tag { display: inline-block; background: #1f6feb33; color: #58a6ff;
         border: 1px solid #1f6feb; border-radius: 12px;
         padding: 2px 10px; font-size: 0.75rem; margin-bottom: 12px; }
  .examples { margin-top: 20px; }
  .examples span { color: #8b949e; font-size: 0.8rem; }
  .ex-btn { background: none; border: 1px solid #30363d; color: #8b949e;
            border-radius: 6px; padding: 4px 10px; margin: 4px;
            cursor: pointer; font-size: 0.8rem; }
  .ex-btn:hover { border-color: #58a6ff; color: #58a6ff; background: none; }
</style>
</head>
<body>
<h1>Mat3ra Prompt-to-Code</h1>
<p>Type a plain-English request and get back the Python API code.</p>
<div class="card">
  <textarea id="prompt" placeholder="e.g. List my Si materials"></textarea>
  <button onclick="ask()">Generate Code</button>

  <div class="examples">
    <span>Try:</span>
    <button class="ex-btn" onclick="fill('List my silicon materials')">List Si materials</button>
    <button class="ex-btn" onclick="fill('Create a new material')">Create material</button>
    <button class="ex-btn" onclick="fill('List my workflows')">List workflows</button>
    <button class="ex-btn" onclick="fill('Submit a new job')">Submit job</button>
    <button class="ex-btn" onclick="fill('Show all my jobs')">List jobs</button>
  </div>

  <div class="label">Detected tool</div>
  <div id="tool" class="tag">—</div>

  <div class="label">Generated Python code</div>
  <pre id="code">—</pre>
</div>

<script>
async function ask() {
  const prompt = document.getElementById('prompt').value.trim();
  if (!prompt) return;
  const r = await fetch('/ask', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({prompt})
  });
  const data = await r.json();
  document.getElementById('tool').textContent = data.tool;
  document.getElementById('code').textContent = data.code;
}
function fill(text) {
  document.getElementById('prompt').value = text;
  ask();
}
document.getElementById('prompt').addEventListener('keydown', e => {
  if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); ask(); }
});
</script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):  # quiet logging
        print(f"  {self.address_string()}  {fmt % args}")

    def _send(self, status: int, content_type: str, body: bytes) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):  # CORS preflight
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.end_headers()

    def do_GET(self):
        if self.path in ("/", "/index.html"):
            self._send(200, "text/html; charset=utf-8", HTML_PLAYGROUND.encode())
        else:
            self._send(404, "text/plain", b"Not found")

    def do_POST(self):
        if self.path != "/ask":
            self._send(404, "text/plain", b"Not found")
            return

        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length)
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            self._send(400, "application/json",
                       json.dumps({"error": "invalid JSON"}).encode())
            return

        prompt = payload.get("prompt", "").strip()
        if not prompt:
            self._send(400, "application/json",
                       json.dumps({"error": "prompt is required"}).encode())
            return

        intent = detect_intent(prompt)
        code   = build_code(intent)

        response = {
            "prompt": prompt,
            "tool":   intent["tool"],
            "args":   {k: v for k, v in intent["args"].items()
                       if k != "owner._id"},
            "code":   code,
        }
        self._send(200, "application/json; charset=utf-8",
                   json.dumps(response, indent=2).encode())


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    HOST, PORT = "localhost", 8888
    server = HTTPServer((HOST, PORT), Handler)
    print(f"\n  Mat3ra Prompt-to-Code server running")
    print(f"  Browser playground : http://{HOST}:{PORT}/")
    print(f"  API endpoint       : POST http://{HOST}:{PORT}/ask")
    print(f"\n  Example:")
    print(f'    curl -s -X POST http://{HOST}:{PORT}/ask \\')
    print(f'         -H "Content-Type: application/json" \\')
    print(f'         -d \'{{"prompt": "List my Si materials"}}\' | python3 -m json.tool')
    print(f"\n  Press Ctrl+C to stop.\n")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n  Server stopped.")
