/* Progressive enhancements for FORD-generated pages. No build step or dependencies. */
(() => {
  "use strict";

  const root = document.documentElement;
  const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)");
  const themeToggle = document.getElementById("cf-theme-toggle");
  const systemTheme = window.matchMedia("(prefers-color-scheme: light)");
  const activeTheme = () => root.dataset.theme ||
    (systemTheme.matches ? "light" : "dark");

  if (themeToggle) {
    const updateThemeLabel = () => {
      const next = activeTheme() === "dark" ? "light" : "dark";
      themeToggle.setAttribute("aria-label", "Switch to " + next + " theme");
      themeToggle.setAttribute("title", "Switch to " + next + " theme");
    };
    themeToggle.addEventListener("click", () => {
      root.dataset.theme = activeTheme() === "dark" ? "light" : "dark";
      try {
        window.localStorage.setItem("constituf-theme", root.dataset.theme);
      } catch (_) {
        /* Local storage may be disabled. */
      }
      updateThemeLabel();
    });
    systemTheme.addEventListener?.("change", updateThemeLabel);
    updateThemeLabel();
  }

  // Indicate the matching primary destination without guessing a parent page.
  const currentPath = window.location.pathname.replace(/\/$/, "/index.html");
  document.querySelectorAll(".cf-links a.nav-link").forEach((link) => {
    const path = new URL(link.href, window.location.href).pathname;
    const isApiPage = /\/(module|type|proc|interface|sourcefile)\//.test(currentPath) &&
      path.endsWith("/lists/modules.html");
    if (path === currentPath || isApiPage) {
      link.setAttribute("aria-current", "page");
    }
  });

  const searchInput = document.getElementById("tipue_search_input");
  if (searchInput) {
    document.addEventListener("keydown", (event) => {
      const target = event.target;
      const editable = target instanceof HTMLElement &&
        (target.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(target.tagName));
      if (event.key !== "/" || event.altKey || event.ctrlKey || event.metaKey || editable) return;
      event.preventDefault();
      const collapsed = document.getElementById("cf-navbar-links");
      if (collapsed && window.innerWidth < 1200 &&
          !collapsed.classList.contains("show") && window.bootstrap?.Collapse) {
        collapsed.addEventListener("shown.bs.collapse", () => searchInput.focus(), {once: true});
        window.bootstrap.Collapse.getOrCreateInstance(collapsed).show();
      } else {
        searchInput.focus();
      }
    });
  }

  // Local outline is built from actual headings; content remains usable without JS.
  const article = document.querySelector(".cf-document-content #text");
  const pageToc = document.getElementById("cf-page-toc");
  const inlineToc = document.getElementById("cf-inline-toc");
  const headings = [];
  const outlineLinks = [];

  if (article) {
    article.querySelectorAll("h2, h3").forEach((heading) => {
      const title = heading.textContent.trim();
      if (!title) return;
      if (!heading.id) {
        const base = title.normalize("NFKD").toLowerCase()
          .replace(/[\u0300-\u036f]/g, "")
          .replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "") || "section";
        let candidate = "cf-" + base;
        let suffix = 2;
        while (document.getElementById(candidate)) candidate = "cf-" + base + "-" + suffix++;
        heading.id = candidate;
      }
      headings.push({element: heading, title});
      const anchor = document.createElement("a");
      anchor.className = "cf-heading-link";
      anchor.href = "#" + encodeURIComponent(heading.id);
      anchor.textContent = "#";
      anchor.setAttribute("aria-label", "Link to " + title);
      heading.appendChild(anchor);
    });

    if (headings.length >= 2 && pageToc && inlineToc) {
      const addOutline = (container) => {
        headings.forEach(({element, title}) => {
          const link = document.createElement("a");
          link.href = "#" + encodeURIComponent(element.id);
          link.textContent = title;
          if (element.tagName === "H3") link.className = "cf-toc-subitem";
          container.appendChild(link);
          outlineLinks.push({link, element});
        });
      };
      addOutline(pageToc.querySelector(".cf-toc-list"));
      addOutline(inlineToc.querySelector(".cf-toc-list"));
      pageToc.hidden = false;
      inlineToc.hidden = false;
      document.querySelector(".cf-document-layout")?.classList.add("cf-has-toc");
    }
  }

  // Filter the atlas by material name, without changing the underlying model data.
  const atlasTiles = [...document.querySelectorAll(".cf-benchmark-atlas .cf-benchmark-tile")];
  if (article && atlasTiles.length >= 12) {
    const controls = document.createElement("div");
    controls.className = "cf-atlas-filter";
    const label = document.createElement("label");
    label.htmlFor = "cf-material-filter";
    label.textContent = "Find a material";
    const input = document.createElement("input");
    input.type = "search";
    input.id = "cf-material-filter";
    input.placeholder = "Filter by material name";
    input.autocomplete = "off";
    input.setAttribute("aria-describedby", "cf-atlas-count");
    const count = document.createElement("span");
    count.id = "cf-atlas-count";
    count.className = "cf-atlas-count";
    count.setAttribute("role", "status");
    count.setAttribute("aria-live", "polite");
    count.textContent = atlasTiles.length + " materials";
    controls.append(label, input, count);
    const resultsHeading = headings.find(({element}) =>
      element.tagName === "H2" && element.textContent.startsWith("Material results"))?.element;
    const insertionPoint = resultsHeading?.nextElementSibling;
    if (insertionPoint) insertionPoint.after(controls);
    else atlasTiles[0].closest(".cf-benchmark-atlas").before(controls);
    input.addEventListener("input", () => {
      const query = input.value.normalize("NFKD").toLocaleLowerCase().trim();
      let matches = 0;
      atlasTiles.forEach((tile) => {
        const name = tile.querySelector("span")?.textContent || tile.textContent;
        const visible = name.normalize("NFKD").toLocaleLowerCase().includes(query);
        tile.hidden = !visible;
        if (visible) matches++;
      });
      count.textContent = matches === atlasTiles.length ?
        atlasTiles.length + " materials" :
        matches === 0 ? "No matching materials" : matches + " of " + atlasTiles.length + " materials";
    });
  }

  // Only show copy controls when the secure Clipboard API is available.
  if (navigator.clipboard?.writeText) {
    document.querySelectorAll(".cf-main pre").forEach((pre) => {
      if (pre.closest(".cf-code-frame")) return;
      const source = pre.textContent;
      if (!source.trim()) return;
      const wrapper = document.createElement("div");
      wrapper.className = "cf-code-frame";
      pre.before(wrapper);
      wrapper.appendChild(pre);
      const button = document.createElement("button");
      button.type = "button";
      button.className = "cf-copy-button";
      button.textContent = "Copy";
      button.setAttribute("aria-label", "Copy code to clipboard");
      button.setAttribute("aria-live", "polite");
      button.addEventListener("click", async () => {
        try {
          await navigator.clipboard.writeText(source);
          button.textContent = "Copied";
          button.setAttribute("aria-label", "Code copied to clipboard");
        } catch (_) {
          button.textContent = "Copy unavailable";
          button.setAttribute("aria-label", "Copy unavailable");
        }
        window.setTimeout(() => {
          button.textContent = "Copy";
          button.setAttribute("aria-label", "Copy code to clipboard");
        }, 1800);
      });
      wrapper.appendChild(button);
    });
  }

  const progress = document.getElementById("cf-reading-progress");
  const backToTop = document.getElementById("cf-back-to-top");
  if (backToTop) {
    backToTop.addEventListener("click", () => {
      window.scrollTo({top: 0, behavior: reducedMotion.matches ? "auto" : "smooth"});
    });
  }
  let pending = false;
  const updateScrollState = () => {
    const range = document.documentElement.scrollHeight - window.innerHeight;
    const fraction = range > 0 ? Math.min(1, Math.max(0, window.scrollY / range)) : 0;
    if (progress) progress.style.transform = "scaleX(" + fraction + ")";
    if (backToTop) backToTop.hidden = window.scrollY < 480;

    if (outlineLinks.length) {
      let active = headings[0].element;
      for (const {element} of headings) {
        if (element.getBoundingClientRect().top <= 148) active = element;
        else break;
      }
      outlineLinks.forEach(({link, element}) => {
        if (element === active) link.setAttribute("aria-current", "location");
        else link.removeAttribute("aria-current");
      });
    }
    pending = false;
  };
  const requestScrollUpdate = () => {
    if (pending) return;
    pending = true;
    window.requestAnimationFrame(updateScrollState);
  };
  window.addEventListener("scroll", requestScrollUpdate, {passive: true});
  window.addEventListener("resize", requestScrollUpdate);
  updateScrollState();
})();
