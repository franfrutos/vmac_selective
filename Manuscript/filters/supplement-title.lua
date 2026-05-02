-- Add a line above the manuscript title on the title page when
-- `supplementary-title-prefix` is defined in metadata.

local stringify = pandoc.utils.stringify

local function text_to_inlines(text)
  local parsed = pandoc.read(text, "markdown")
  if parsed.blocks[1] and parsed.blocks[1].content then
    return parsed.blocks[1].content
  end
  return pandoc.List:new()
end

function Meta(meta)
  if not meta["supplementary-title-prefix"] or not meta.title then
    return meta
  end

  local prefix = text_to_inlines(stringify(meta["supplementary-title-prefix"]))
  local title = text_to_inlines(stringify(meta.title))

  if FORMAT:match("latex") then
    local combined = pandoc.List:new()
    combined:extend(prefix)
    combined:insert(pandoc.RawInline("latex", "\\\\[1em]"))
    combined:extend(title)
    meta.title = pandoc.MetaInlines(combined)
    return meta
  end

  local combined = pandoc.List:new()
  combined:extend(prefix)
  combined:insert(pandoc.LineBreak())
  combined:insert(pandoc.LineBreak())
  combined:extend(title)
  meta.title = pandoc.MetaInlines(combined)
  return meta
end
