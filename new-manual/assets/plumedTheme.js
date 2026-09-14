(function () {
  "use strict";

  var root = document.documentElement;
  var themes = ["auto", "light", "dark"];
  var savedTheme = localStorage.getItem("plumed-doc-theme") || "auto";
  if (themes.indexOf(savedTheme) === -1) savedTheme = "auto";
  root.dataset.theme = savedTheme;

  var themeButton = document.getElementById("theme-toggle");
  themeButton.addEventListener("click", function () {
    var next = themes[(themes.indexOf(root.dataset.theme) + 1) % themes.length];
    root.dataset.theme = next;
    localStorage.setItem("plumed-doc-theme", next);
    themeButton.textContent = "Theme: " + next;
  });
  themeButton.textContent = "Theme: " + savedTheme;

  var menuButton = document.querySelector(".menu-toggle");
  var navigation = document.getElementById("site-nav");
  menuButton.addEventListener("click", function () {
    var open = navigation.classList.toggle("open");
    menuButton.setAttribute("aria-expanded", String(open));
  });

  document.querySelectorAll("blockquote > p:first-child").forEach(function (title) {
    var match = title.textContent.match(/^\[!([A-Z_-]+)\]\s*(.*)$/);
    if (!match) return;
    var block = title.parentElement;
    block.classList.add("admonition", match[1].toLowerCase());
    title.classList.add("admonition-title");
    title.textContent = match[2] || match[1].toLowerCase().replace(/_/g, " ");
  });

  var dialog = document.getElementById("search-dialog");
  var query = document.getElementById("search-query");
  var results = document.getElementById("search-results");
  var searchIndex = null;

  document.getElementById("search-open").addEventListener("click", function () {
    dialog.showModal();
    query.focus();
    if (!searchIndex) {
      fetch(document.body.dataset.searchIndex)
        .then(function (response) { return response.json(); })
        .then(function (index) { searchIndex = index; runSearch(); });
    }
  });

  function runSearch() {
    if (!searchIndex) return;
    var terms = query.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
    results.replaceChildren();
    if (!terms.length) return;
    searchIndex.filter(function (entry) {
      var haystack = (entry.title + " " + entry.text).toLowerCase();
      return terms.every(function (term) { return haystack.indexOf(term) !== -1; });
    }).slice(0, 30).forEach(function (entry) {
      var item = document.createElement("li");
      var link = document.createElement("a");
      var rootPath = document.body.dataset.siteRoot.replace(/\/$/, "");
      link.href = rootPath + entry.location;
      link.textContent = entry.title;
      var excerpt = document.createElement("small");
      excerpt.textContent = entry.text.slice(0, 180) + (entry.text.length > 180 ? "…" : "");
      item.append(link, excerpt);
      results.append(item);
    });
  }

  query.addEventListener("input", runSearch);
}());
