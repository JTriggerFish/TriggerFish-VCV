#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace tfseq::editor {

struct EditResult {
  std::string text;
  int selection = 0;
  int cursor = 0;
  bool changed = false;
};

inline int ClampPosition(const std::string &text, int position) noexcept {
  return std::clamp(position, 0, static_cast<int>(text.size()));
}

inline bool IsLineBreak(const char character) noexcept {
  return character == '\n' || character == '\r';
}

inline int LineBegin(const std::string &text, int position) noexcept {
  position = ClampPosition(text, position);
  if (position == 0)
    return 0;
  int search = position - 1;
  // A cursor on the LF half of CRLF still belongs to the preceding line.
  if (position < static_cast<int>(text.size()) &&
      text[static_cast<std::size_t>(position)] == '\n' && search >= 0 &&
      text[static_cast<std::size_t>(search)] == '\r')
    --search;
  for (int index = search; index >= 0; --index) {
    const char character = text[static_cast<std::size_t>(index)];
    if (character == '\n')
      return index + 1;
    if (character == '\r')
      return index +
             (index + 1 < static_cast<int>(text.size()) &&
                      text[static_cast<std::size_t>(index + 1)] == '\n'
                  ? 2
                  : 1);
  }
  return 0;
}

inline int LineEnd(const std::string &text, int position) noexcept {
  position = ClampPosition(text, position);
  int search = position;
  if (search < static_cast<int>(text.size()) &&
      text[static_cast<std::size_t>(search)] == '\n' && search > 0 &&
      text[static_cast<std::size_t>(search - 1)] == '\r')
    --search;
  for (int index = search; index < static_cast<int>(text.size()); ++index) {
    if (IsLineBreak(text[static_cast<std::size_t>(index)]))
      return index;
  }
  return static_cast<int>(text.size());
}

inline int NextLineBegin(const std::string &text, int lineEnd) noexcept {
  lineEnd = ClampPosition(text, lineEnd);
  if (lineEnd >= static_cast<int>(text.size()))
    return lineEnd;
  if (text[static_cast<std::size_t>(lineEnd)] == '\r' &&
      lineEnd + 1 < static_cast<int>(text.size()) &&
      text[static_cast<std::size_t>(lineEnd + 1)] == '\n')
    return lineEnd + 2;
  return lineEnd + 1;
}

inline std::string LineBreakAt(const std::string &text, int lineEnd) {
  const int next = NextLineBegin(text, lineEnd);
  if (next > lineEnd)
    return text.substr(static_cast<std::size_t>(lineEnd),
                       static_cast<std::size_t>(next - lineEnd));
  return {};
}

inline std::string PreferredLineBreak(const std::string &text,
                                      int lineBegin) {
  if (lineBegin >= 2 && text[static_cast<std::size_t>(lineBegin - 2)] == '\r' &&
      text[static_cast<std::size_t>(lineBegin - 1)] == '\n')
    return "\r\n";
  if (lineBegin >= 1 &&
      IsLineBreak(text[static_cast<std::size_t>(lineBegin - 1)]))
    return std::string(1, text[static_cast<std::size_t>(lineBegin - 1)]);
  for (int index = 0; index < static_cast<int>(text.size()); ++index) {
    if (!IsLineBreak(text[static_cast<std::size_t>(index)]))
      continue;
    if (text[static_cast<std::size_t>(index)] == '\r' &&
        index + 1 < static_cast<int>(text.size()) &&
        text[static_cast<std::size_t>(index + 1)] == '\n')
      return "\r\n";
    return std::string(1, text[static_cast<std::size_t>(index)]);
  }
  return "\n";
}

inline EditResult ToggleLineComments(const std::string &source, int cursor,
                                     int selection) {
  EditResult result{source, selection, cursor, false};
  cursor = ClampPosition(source, cursor);
  selection = ClampPosition(source, selection);
  const int selectedBegin = std::min(cursor, selection);
  const int selectedEnd = std::max(cursor, selection);
  const int blockBegin = LineBegin(source, selectedBegin);
  int finalSelectedPosition = selectedEnd;
  if (selectedEnd > selectedBegin && selectedEnd > 0 &&
      IsLineBreak(source[static_cast<std::size_t>(selectedEnd - 1)]))
    --finalSelectedPosition;
  const int blockEnd = LineEnd(source, finalSelectedPosition);

  struct LineEdit {
    int position = 0;
    int remove = 0;
    std::string insert;
  };
  std::vector<LineEdit> edits;
  bool allCommented = true;
  bool hasContent = false;
  for (int line = blockBegin; line <= blockEnd;) {
    const int end = LineEnd(source, line);
    int content = line;
    while (content < end && (source[static_cast<std::size_t>(content)] == ' ' ||
                             source[static_cast<std::size_t>(content)] == '\t'))
      ++content;
    if (content < end) {
      hasContent = true;
      if (content + 1 >= end ||
          source.compare(static_cast<std::size_t>(content), 2, "//") != 0)
        allCommented = false;
    }
    if (end >= static_cast<int>(source.size()))
      break;
    line = NextLineBegin(source, end);
  }
  if (!hasContent)
    return result;

  for (int line = blockBegin; line <= blockEnd;) {
    const int end = LineEnd(source, line);
    int content = line;
    while (content < end && (source[static_cast<std::size_t>(content)] == ' ' ||
                             source[static_cast<std::size_t>(content)] == '\t'))
      ++content;
    if (content < end) {
      if (allCommented) {
        int remove = 2;
        if (content + remove < end &&
            source[static_cast<std::size_t>(content + remove)] == ' ')
          ++remove;
        edits.push_back({content, remove, {}});
      } else {
        edits.push_back({content, 0, "// "});
      }
    }
    if (end >= static_cast<int>(source.size()))
      break;
    line = NextLineBegin(source, end);
  }

  int delta = 0;
  for (const auto &edit : edits)
    delta += static_cast<int>(edit.insert.size()) - edit.remove;
  for (auto edit = edits.rbegin(); edit != edits.rend(); ++edit)
    result.text.replace(static_cast<std::size_t>(edit->position),
                        static_cast<std::size_t>(edit->remove), edit->insert);
  result.selection = blockBegin;
  result.cursor = blockEnd + delta;
  result.changed = true;
  return result;
}

inline EditResult Duplicate(const std::string &source, int cursor,
                            int selection) {
  EditResult result{source, selection, cursor, false};
  cursor = ClampPosition(source, cursor);
  selection = ClampPosition(source, selection);
  if (cursor != selection) {
    const int begin = std::min(cursor, selection);
    const int end = std::max(cursor, selection);
    const std::string copy = source.substr(
        static_cast<std::size_t>(begin), static_cast<std::size_t>(end - begin));
    result.text.insert(static_cast<std::size_t>(end), copy);
    result.selection = end;
    result.cursor = end + static_cast<int>(copy.size());
    result.changed = true;
    return result;
  }
  if (source.empty())
    return result;

  const int begin = LineBegin(source, cursor);
  const int end = LineEnd(source, cursor);
  const int column = std::clamp(cursor - begin, 0, end - begin);
  const std::string line = source.substr(static_cast<std::size_t>(begin),
                                         static_cast<std::size_t>(end - begin));
  if (end < static_cast<int>(source.size())) {
    const int insertion = NextLineBegin(source, end);
    result.text.insert(static_cast<std::size_t>(insertion),
                       line + LineBreakAt(source, end));
    result.cursor = insertion + column;
  } else {
    const std::string lineBreak = PreferredLineBreak(source, begin);
    result.text.insert(static_cast<std::size_t>(end), lineBreak + line);
    result.cursor = end + static_cast<int>(lineBreak.size()) + column;
  }
  result.selection = result.cursor;
  result.changed = true;
  return result;
}

} // namespace tfseq::editor
