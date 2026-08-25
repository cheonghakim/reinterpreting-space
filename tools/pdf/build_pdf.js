// Renders lapse_gradient_topology_manuscript_en.md to a PDF with real
// KaTeX-rendered math (the previous PDF export pipeline printed raw "$$
// \lambda ... $$" text because it never invoked a math renderer).
//
// Usage (from this directory, after `npm install`):
//   node build_pdf.js [source.md] [output.pdf]
// Both arguments are optional; they default to the manuscript at the repo
// root. CHROME_PATH can override browser auto-detection.

const fs = require('fs');
const path = require('path');
const { pathToFileURL } = require('url');
const os = require('os');
const MarkdownIt = require('markdown-it');
const texmath = require('markdown-it-texmath');
const katex = require('katex');
const puppeteer = require('puppeteer-core');

const REPO_ROOT = path.resolve(__dirname, '..', '..');
const SRC = path.resolve(process.argv[2] || path.join(REPO_ROOT, 'lapse_gradient_topology_manuscript_en.md'));
const OUT_PDF = path.resolve(process.argv[3] || SRC.replace(/\.md$/, '.pdf'));
const DOCDIR = path.dirname(SRC);
const KATEX_CSS = path.join(path.dirname(require.resolve('katex/package.json')), 'dist', 'katex.min.css');

function findBrowser() {
  if (process.env.CHROME_PATH && fs.existsSync(process.env.CHROME_PATH)) return process.env.CHROME_PATH;
  const candidates = [
    'C:/Program Files/Google/Chrome/Application/chrome.exe',
    'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe',
    'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
    'C:/Program Files/Microsoft/Edge/Application/msedge.exe',
    '/usr/bin/google-chrome',
    '/usr/bin/chromium-browser',
    '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  ];
  const found = candidates.find((p) => fs.existsSync(p));
  if (!found) {
    console.error(
      'No local Chrome/Edge found in the usual install locations.\n' +
      'Set CHROME_PATH to your browser executable and re-run, e.g.\n' +
      '  CHROME_PATH="C:/Program Files/Google/Chrome/Application/chrome.exe" node build_pdf.js'
    );
    process.exit(1);
  }
  return found;
}

function render() {
  const md = new MarkdownIt({ html: true, linkify: true, typographer: true });
  md.use(texmath, {
    engine: katex,
    delimiters: ['dollars', 'brackets'], // supports both $$...$$ and \( \) \[ \]
    katexOptions: { throwOnError: false, macros: {} },
  });

  const src = fs.readFileSync(SRC, 'utf8');
  const body = md.render(src);

  return `<!doctype html>
<html>
<head>
<meta charset="utf-8">
<base href="${pathToFileURL(DOCDIR + path.sep).href}">
<link rel="stylesheet" href="${pathToFileURL(KATEX_CSS).href}">
<style>
  body {
    font-family: "Georgia", "Times New Roman", serif;
    font-size: 11.5pt;
    line-height: 1.55;
    color: #111;
    max-width: 720px;
    margin: 0 auto;
    padding: 0 8px;
  }
  h1 { font-size: 20pt; margin-top: 0; }
  h2 { font-size: 15pt; border-bottom: 1px solid #ccc; padding-bottom: 4px; margin-top: 2em; }
  h3 { font-size: 12.5pt; margin-top: 1.6em; }
  code, pre { font-family: "Consolas", "Menlo", monospace; font-size: 0.92em; }
  pre { background: #f5f5f5; padding: 10px; border-radius: 4px; overflow-x: auto; }
  code { background: #f0f0f0; padding: 1px 4px; border-radius: 3px; }
  pre code { background: none; padding: 0; }
  blockquote { border-left: 3px solid #999; margin-left: 0; padding-left: 1em; color: #444; }
  img { max-width: 100%; display: block; margin: 1em auto; }
  hr { border: none; border-top: 1px solid #ccc; margin: 2em 0; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; }
  th, td { border: 1px solid #ccc; padding: 6px 10px; text-align: left; }
  .katex-display { margin: 0.9em 0; overflow-x: auto; overflow-y: hidden; }
  .katex { font-size: 1.05em; }
</style>
</head>
<body>
${body}
</body>
</html>`;
}

(async () => {
  if (!fs.existsSync(SRC)) {
    console.error('Source Markdown file not found:', SRC);
    process.exit(1);
  }
  const chromePath = findBrowser();
  const html = render();
  const tmpHtml = path.join(os.tmpdir(), 'validated_mp_manuscript_render.html');
  fs.writeFileSync(tmpHtml, html, 'utf8');

  const browser = await puppeteer.launch({ executablePath: chromePath, headless: 'new' });
  try {
    const page = await browser.newPage();
    await page.goto(pathToFileURL(tmpHtml).href, { waitUntil: 'networkidle0' });
    await page.pdf({
      path: OUT_PDF,
      format: 'A4',
      margin: { top: '20mm', bottom: '20mm', left: '18mm', right: '18mm' },
      printBackground: true,
    });
  } finally {
    await browser.close();
  }
  console.log('Wrote', OUT_PDF);
})().catch((err) => {
  console.error(err);
  process.exit(1);
});
