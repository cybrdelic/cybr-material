// Run against a local repository-root HTTP server. Playwright is optional tooling.
// GALLERY_BROWSER=/path/to/chrome node tools/gallery/browser-check.cjs
const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');

(async () => {
  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.GALLERY_BROWSER || undefined,
    args: ['--no-sandbox', '--disable-dev-shm-usage'],
  });
  const base = process.env.GALLERY_ORIGIN || 'http://127.0.0.1:8000';
  const screenshots = process.env.GALLERY_SCREENSHOTS;
  const report = [];
  try {
    for (const viewport of [{ width: 390, height: 844 }, { width: 1440, height: 1000 }]) {
      const context = await browser.newContext({ viewport, deviceScaleFactor: 1 });
      const page = await context.newPage();
      const errors = [];
      const requests = [];
      let receivedBytes = 0;
      page.on('pageerror', error => errors.push(error.message));
      page.on('request', request => requests.push(request.url()));
      page.on('response', response => {
        if (response.status() >= 400) errors.push(`${response.status()} ${response.url()}`);
      });
      const cdp = await context.newCDPSession(page);
      await cdp.send('Network.enable');
      cdp.on('Network.loadingFinished', event => { receivedBytes += event.encodedDataLength; });
      await page.goto(`${base}/gallery/`, { waitUntil: 'networkidle' });
      await page.locator('#material-grid .card img').first().waitFor();
      const initialBytes = receivedBytes;
      assert(initialBytes < 1024 * 1024, `Initial document exceeded 1 MiB: ${initialBytes}`);
      assert(!requests.some(url => url.endsWith('.png')), 'Native PNG requested before explicit navigation');
      assert.equal(await page.locator('#material-grid .card:visible').count(), 10);
      await page.selectOption('#family', 'Wood');
      assert.equal(await page.locator('#material-grid .card:visible').count(), 2);
      assert.equal(await page.locator('#result-count').innerText(), '2 materials');
      await page.selectOption('#family', 'all');
      assert.equal(await page.locator('#material-grid .card:visible').count(), 10);
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Home overflow');
      if (screenshots) {
        await fs.mkdir(screenshots, { recursive: true });
        await page.locator('img').evaluateAll(images => images.forEach(image => { image.loading = 'eager'; }));
        await page.waitForLoadState('networkidle');
        await page.screenshot({ path: path.join(screenshots, `gallery-${viewport.width}.png`), fullPage: true });
      }
      await page.locator('#material-grid a').first().click();
      await page.waitForLoadState('networkidle');
      assert.equal(await page.locator('.channel').count(), 12);
      assert.equal(await page.locator('.view').count(), 2);
      assert(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'Material page overflow');
      const nativeLink = page.locator('.channel').first().getByRole('link', { name: 'Open native 4K PNG' });
      const href = await nativeLink.getAttribute('href');
      assert(href.endsWith('/BaseColor.png'));
      const response = await page.request.get(new URL(href, page.url()).href);
      assert(response.ok());
      const bytes = await response.body();
      assert.equal(bytes.readUInt32BE(16), 4096);
      assert.equal(bytes.readUInt32BE(20), 4096);
      if (screenshots) {
        await page.locator('img').evaluateAll(images => images.forEach(image => { image.loading = 'eager'; }));
        await page.waitForLoadState('networkidle');
        await page.screenshot({ path: path.join(screenshots, `material-${viewport.width}.png`), fullPage: true });
      }
      assert.deepEqual(errors, []);
      report.push({ viewport, initialBytes, nativePngSize: [4096, 4096], overflow: false, errors });
      await context.close();
    }
    const noJs = await browser.newContext({ javaScriptEnabled: false, viewport: { width: 390, height: 844 } });
    const page = await noJs.newPage();
    await page.goto(`${base}/gallery/`);
    assert.equal(await page.locator('#material-grid a').count(), 10);
    assert.equal(await page.locator('.filters:visible').count(), 0);
    await page.locator('#material-grid a').first().click();
    assert.equal(await page.locator('.channel').count(), 12);
    await noJs.close();
    console.log(JSON.stringify({ result: 'PASS', report, withoutJavaScript: 'native navigation and channels available' }, null, 2));
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
