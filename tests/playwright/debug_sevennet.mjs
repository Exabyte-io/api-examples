import { chromium } from "playwright";
const URL = "http://localhost:8000/lab/index.html?path=experiments/relax_structure_with_sevennet.ipynb";
const b = await chromium.launch({ channel: "chrome", headless: false });
const p = await b.newPage();
await p.goto(URL, { waitUntil: "domcontentloaded", timeout: 60000 });
await p.waitForFunction(() => {
  const ind = document.querySelector(".jp-Notebook-ExecutionIndicator");
  return ind?.dataset?.status === "idle" && document.querySelectorAll(".jp-CodeCell").length > 5;
}, null, { timeout: 120000, polling: 2000 });
await p.locator('button[data-command="runmenu:restart-and-run-all"]').click();
await p.waitForTimeout(1000);
try { await p.locator(".jp-Dialog button", { hasText: /restart/i }).click({ timeout: 3000 }); } catch {}
for (let i = 0; i < 30; i++) {
  await p.waitForTimeout(5000);
  const s = await p.evaluate(() => {
    const k = document.querySelector(".jp-Notebook-ExecutionIndicator")?.dataset?.status;
    let running = 0;
    document.querySelectorAll(".jp-CodeCell .jp-InputPrompt").forEach(e => {
      if (e.textContent.trim() === "[*]:") running++;
    });
    return { k, running };
  });
  if (s.k === "idle" && s.running === 0) break;
}
const err = await p.evaluate(() => {
  const cells = document.querySelectorAll(".jp-CodeCell");
  const cell = cells[6]; // cell 7
  if (!cell) return "no cell 7";
  // Get just stderr/error outputs, skip plotly JS
  let text = "";
  cell.querySelectorAll(".jp-OutputArea-output").forEach((e, i) => {
    const t = e.textContent || "";
    // Skip massive plotly JS, only get error outputs
    if (t.includes("Traceback") || t.includes("Error") || t.includes("Step:")) {
      // Get last 3000 chars to capture the actual error
      text += `[output ${i}]: ${t.slice(-3000)}\n`;
    }
  });
  return text || "No error outputs found";
});
console.log(err);
await b.close();
