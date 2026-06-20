
const toggle = document.querySelector('.nav-toggle');
const nav = document.querySelector('.site-nav');
if (toggle && nav) {
  toggle.addEventListener('click', () => {
    const open = nav.classList.toggle('open');
    toggle.setAttribute('aria-expanded', String(open));
  });
}

// Ensure external and document links open in a new tab.
document.querySelectorAll('a[href]').forEach((link) => {
  const href = link.getAttribute('href') || '';
  const isDocument = /\.(pdf|docx?|pptx?|xlsx?)$/i.test(href);
  const isExternal = href.startsWith('http://') || href.startsWith('https://');
  if (isDocument || isExternal) {
    link.setAttribute('target', '_blank');
    link.setAttribute('rel', 'noopener noreferrer');
  }
});
