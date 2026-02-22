--- custom-keywords.lua
--- Reads `custom-keywords` from YAML metadata and appends them after the
--- abstract block produced by the abstract-section filter.
--- Each entry has `name` (e.g. "JEL Codes") and `values` (a list).
--- For LaTeX output the hikmah `keywords` environment style is reused;
--- for HTML a simple paragraph is emitted.

local function render_latex(groups)
  local lines = {}
  for _, grp in ipairs(groups) do
    local name = pandoc.utils.stringify(grp.name)
    local vals = {}
    for _, v in ipairs(grp.values) do
      table.insert(vals, pandoc.utils.stringify(v))
    end
    local joined = table.concat(vals, "; ")
    table.insert(lines,
      "\\vskip -3em \\hspace{\\parindent}" ..
      "{\\sffamily\\footnotesize\\bfseries\\MakeUppercase{" .. name .. "}}\\quad " ..
      "{\\sffamily\\small " .. joined .. "}" ..
      "\\vskip 3em")
  end
  return pandoc.RawBlock("latex", table.concat(lines, "\n"))
end

local function render_html(groups)
  local parts = {}
  for _, grp in ipairs(groups) do
    local name = pandoc.utils.stringify(grp.name)
    local vals = {}
    for _, v in ipairs(grp.values) do
      table.insert(vals, pandoc.utils.stringify(v))
    end
    table.insert(parts,
      "<p><strong>" .. name .. ":</strong> " ..
      table.concat(vals, "; ") .. "</p>")
  end
  return pandoc.RawBlock("html", table.concat(parts, "\n"))
end

function Pandoc(doc)
  local ck = doc.meta["custom-keywords"]
  if not ck then return nil end

  -- Parse the metadata list
  local groups = {}
  for _, item in ipairs(ck) do
    local name = item.name
    local values = item.values or {}
    table.insert(groups, { name = name, values = values })
  end
  if #groups == 0 then return nil end

  -- Determine output format
  local is_latex = FORMAT:match("latex") or FORMAT:match("pdf")
  local block
  if is_latex then
    block = render_latex(groups)
  else
    block = render_html(groups)
  end

  -- Find the abstract header and insert after abstract content.
  -- The abstract-section filter converts the Abstract heading into metadata,
  -- so in the final document the abstract appears via the template.
  -- We insert our block at the very beginning of the body so it appears
  -- right after the template-rendered abstract + keywords.
  table.insert(doc.blocks, 1, block)

  return doc
end
