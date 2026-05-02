-- figletter_link.lua (no-space version)
-- Merge [Fig. 2](#fig-2) + [A]{.figletter} -> "Fig. 2A" (same target), WITHOUT a space.

-- Helpers -------------------------------------------------------------

local function has_class(el, class)
  if not el or el.t ~= "Span" or not el.classes then return false end
  for _, c in ipairs(el.classes) do
    if c == class then return true end
  end
  return false
end

-- Trim trailing whitespace tokens (Space/SoftBreak/LineBreak) and
-- strip trailing spaces (incl. NBSP) inside last Str
local function trim_trailing_ws(inls)
  while #inls > 0 do
    local last = inls[#inls]
    if last.t == "Space" or last.t == "SoftBreak" or last.t == "LineBreak" then
      table.remove(inls) -- pop
    elseif last.t == "Str" then
      -- remove trailing ASCII/Unicode NBSP
      local txt = last.text
      -- strip spaces and NBSP (\u{00A0})
      txt = txt:gsub("[ %z\194\160]+$", "")
      if txt == "" then
        table.remove(inls)
      else
        last.text = txt
        break
      end
    else
      break
    end
  end
end

-- Trim leading whitespace tokens and strip leading spaces in first Str
local function trim_leading_ws(inls)
  while #inls > 0 do
    local first = inls[1]
    if first.t == "Space" or first.t == "SoftBreak" or first.t == "LineBreak" then
      table.remove(inls, 1)
    elseif first.t == "Str" then
      local txt = first.text
      txt = txt:gsub("^[ %z\194\160]+", "")
      if txt == "" then
        table.remove(inls, 1)
      else
        first.text = txt
        break
      end
    else
      break
    end
  end
end

return {
  Inlines = function(inlines)
    local i = 1
    while i <= #inlines do
      local link = inlines[i]

      if link and link.t == "Link" and type(link.target) == "string"
         and link.target:match("^#fig[%-%:]") then

        local j = i + 1

        -- Remove ANY number of Space/SoftBreak directly after the link
        while inlines[j] and (inlines[j].t == "Space" or inlines[j].t == "SoftBreak") do
          table.remove(inlines, j)
        end

        local span = inlines[j]
        if has_class(span, "figletter") then
          -- Ensure no trailing/leading whitespace at the join
          trim_trailing_ws(link.content)
          trim_leading_ws(span.content)

          -- Append span content (e.g., "A") right after the link text
          for _, el in ipairs(span.content) do
            table.insert(link.content, el)
          end

          -- Remove the span
          table.remove(inlines, j)
          -- Nothing more to remove; any gap was already stripped above
        end
      end

      i = i + 1
    end
    return inlines
  end
}
