-- Apply blue highlighting for selected inline text in PDF/LaTeX and DOCX.
--
-- Usage in markdown/quarto:
--   [this text will be blue]{.blue}
--   [this also works]{.review-blue}
--   [and this too]{.azul}

local TARGET_CLASSES = {
  blue = true,
  ["review-blue"] = true,
  azul = true
}

local function has_target_class(el)
  if not el or not el.classes then
    return false
  end

  for _, class in ipairs(el.classes) do
    if TARGET_CLASSES[class] then
      return true
    end
  end

  return false
end

local function ensure_xcolor(meta)
  if not FORMAT:match("latex") then
    return meta
  end

  local include = meta["include-in-header"]
  local package = pandoc.RawBlock("latex", "\\usepackage{xcolor}")

  if not include then
    meta["include-in-header"] = pandoc.MetaBlocks({ package })
    return meta
  end

  if include.t == "MetaBlocks" then
    include[#include + 1] = package
    meta["include-in-header"] = include
    return meta
  end

  if include.t == "MetaInlines" then
    meta["include-in-header"] = pandoc.MetaBlocks({
      pandoc.Plain(include),
      package
    })
    return meta
  end

  if include.t == "MetaList" then
    include[#include + 1] = pandoc.MetaBlocks({ package })
    meta["include-in-header"] = include
    return meta
  end

  return meta
end

function Meta(meta)
  return ensure_xcolor(meta)
end

function Span(el)
  if not has_target_class(el) then
    return nil
  end

  if FORMAT == "docx" then
    el.attributes["custom-style"] = "Hyperlink"
    return el
  end

  if FORMAT:match("latex") then
    local colored = { pandoc.RawInline("latex", "\\textcolor{blue}{") }
    for _, inline in ipairs(el.content) do
      colored[#colored + 1] = inline
    end
    colored[#colored + 1] = pandoc.RawInline("latex", "}")
    return colored
  end

  return el
end
