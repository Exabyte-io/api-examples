# Mat3ra MCP Server

Proof-of-concept [Model Context Protocol](https://modelcontextprotocol.io/) server
that exposes the Mat3ra REST API as MCP **tools**, allowing any compatible client
(Claude Desktop, Cursor, etc.) to interact with your Mat3ra account through natural language.

## Available tools

| Tool | Description |
|---|---|
| `list_materials` | Search materials by formula or Mongo query |
| `get_material` | Fetch one material by ID |
| `create_material` | Upload a new material from a JSON config |
| `list_workflows` | List workflows |
| `list_jobs` | List jobs |
| `create_and_submit_job` | Create + submit a job (requires material + workflow IDs) |
| `get_job` | Fetch a job and its current status |

## Setup

```bash
# 1. Create a venv (Python 3.10+)
python3.11 -m venv .venv
source .venv/bin/activate

# 2. Install dependencies
pip install mcp mat3ra-api-client

# 3. Set your credentials
cp mcp_server/.env.example mcp_server/.env
#   → edit .env with your ACCOUNT_ID and AUTH_TOKEN
#   (get them from https://platform.mat3ra.com → Account Preferences → API Tokens)

# 4. Start the server (stdio transport)
export $(cat mcp_server/.env | grep -v '^#' | xargs)
python mcp_server/server.py
```

## Claude Desktop integration

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "mat3ra": {
      "command": "/path/to/.venv/bin/python",
      "args": ["/path/to/api-examples/mcp_server/server.py"],
      "env": {
        "MAT3RA_ACCOUNT_ID": "your_account_id",
        "MAT3RA_AUTH_TOKEN": "your_auth_token"
      }
    }
  }
}
```

Then restart Claude Desktop.  You can now ask things like:
- *"List my Si materials on Mat3ra"*
- *"Create an FCC Si material with lattice parameter 5.43 Å"*
- *"What is the status of job `abc123`?"*
