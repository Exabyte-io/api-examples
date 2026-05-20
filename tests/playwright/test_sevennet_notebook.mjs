import { chromium } from "playwright";

const URL =
  "http://localhost:8000/lab/index.html?path=experiments/relax_structure_with_sevennet.ipynb";

(async () => {
  const browser = await chromium.launch({ channel: "chrome", headless: false });
  const page = await browser.newPage();

  console.log("⏳ Navigating to JupyterLite...");
  await page.goto(URL, { waitUntil: "domcontentloaded", timeout: 60_000 });

  // Wait for kernel idle
  console.log("⏳ Waiting for Pyodide kernel...");
  await page.waitForFunction(
    () => {
      const ind = document.querySelector(".jp-Notebook-ExecutionIndicator");
      return ind?.dataset?.status === "idle" && document.querySelectorAll(".jp-CodeCell").length > 5;
    },
    null,
    { timeout: 120_000, polling: 2_000 }
  );
  console.log("✅ Kernel idle!");

  // Click Restart & Run All
  await page.locator('button[data-command="runmenu:restart-and-run-all"]').click();
  await page.waitForTimeout(1_000);
  try {
    await page.locator(".jp-Dialog button", { hasText: /restart/i }).click({ timeout: 3_000 });
    console.log("✅ Restart & Run All confirmed!");
  } catch {
    console.log("⚠ No dialog");
  }

  // Monitor with logging
  console.log("\n📊 Monitoring...");
  const start = Date.now();

  while (Date.now() - start < 300_000) {
    await page.waitForTimeout(10_000);

    const s = await page.evaluate(() => {
      const cells = document.querySelectorAll(".jp-CodeCell");
      let running = 0, done = 0, errors = 0;
      cells.forEach((c) => {
        const p = c.querySelector(".jp-InputPrompt")?.textContent?.trim() || "";
        if (p === "[*]:") running++;
        else if (/\[\d+\]:/.test(p)) done++;
        if (c.querySelector(".jp-RenderedText[data-mime-type='application/vnd.jupyter.stderr']")) errors++;
      });
      const k = document.querySelector(".jp-Notebook-ExecutionIndicator")?.dataset?.status;
      return { running, done, errors, total: cells.length, k };
    });

    const t = ((Date.now() - start) / 1000) | 0;
    console.log(`  [${t}s] kernel=${s.k} running=${s.running} done=${s.done}/${s.total} errors=${s.errors}`);

    if (s.k === "idle" && s.running === 0 && s.done > 0) break;
  }

  // Final report
  console.log("\n" + "=".repeat(60));
  const results = await page.evaluate(() =>
    Array.from(document.querySelectorAll(".jp-CodeCell")).map((c, i) => {
      let t = "";
      for (const o of c.querySelectorAll(".jp-OutputArea-output")) t += o.textContent;
      return { cell: i + 1, err: t.includes("Traceback"), out: t };
    })
  );

  let failures = 0;
  for (const r of results) {
    if (r.err) {
      failures++;
      console.log(`Cell ${r.cell} ❌`);
      console.log(r.out.substring(0, 2000));
      console.log();
    } else {
      const preview = r.out.trim().substring(0, 120).replace(/\n/g, " | ");
      console.log(`Cell ${r.cell} ✅${preview ? ": " + preview : ""}`);
    }
  }

  console.log("=".repeat(60));
  if (failures > 0) {
    console.log(`\n⚠ ${failures} cell(s) failed`);
    process.exitCode = 1;
  } else {
    console.log(`\n🎉 ALL ${results.length} CELLS PASSED!`);
  }

  await browser.close();
})();
