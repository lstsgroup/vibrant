from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.statemachine import StringList


class KeywordNode(nodes.General, nodes.Element):
    """Custom node for input keyword entries."""
    pass


def visit_keyword_node_html(self, node):
    section = node.get('section', '')
    keyword = node.get('keyword', '')
    self.body.append(
        f"<div class='keyword-entry'>"
        f"<div class='keyword-header'>{section} :: <strong>{keyword}</strong></div>"
    )


def depart_keyword_node_html(self, node):
    self.body.append("</div>")  # closes .keyword-entry


class KeywordDirective(Directive):
    required_arguments = 1
    option_spec = {
        'section': directives.unchanged_required,
        'type': directives.unchanged,
        'default': directives.unchanged,
        'unit': directives.unchanged,
    }
    has_content = True

    def run(self):
        keyword = self.arguments[0]
        section = self.options.get('section', 'general')
        type_ = self.options.get('type', '—')
        default = self.options.get('default', '—')
        unit = self.options.get('unit', None)
        description = '\n'.join(self.content)

        # Create a target for cross-referencing
        target_id = f'keyword-{keyword.lower()}'
        target_node = nodes.target('', '', ids=[target_id])

        # Store info in the node (for HTML visitor)
        node = KeywordNode()
        node['section'] = section
        node['keyword'] = keyword

        # Build reST for type/default/unit/description
        text = [
            f"- **Type**: {type_}",
            f"- **Default**: {default}",
        ]
        if unit:
            text.append(f"- **Unit**: {unit}")
        text += ["-  \u200b", "- " + description]

        content = StringList(text)
        self.state.nested_parse(content, self.content_offset, node)

        return [target_node, node]


def setup(app):
    app.add_node(
        KeywordNode,
        html=(visit_keyword_node_html, depart_keyword_node_html)
    )
    app.add_directive("keyword", KeywordDirective)
