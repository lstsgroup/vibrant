from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.statemachine import StringList

class KeywordNode(nodes.General, nodes.Element):
    """Custom node for input keyword entries."""
    pass


def visit_keyword_node_html(self, node):
    self.body.append(self.starttag(node, 'div', CLASS='keyword-entry'))


def depart_keyword_node_html(self, node):
    self.body.append('</div>')


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

        # Build a unique ID for referencing
        target_id = f'keyword-{keyword.lower()}'
        target_node = nodes.target('', '', ids=[target_id])

        # Build text block for nested parsing
        text = [
            f"**{section} :: {keyword}**",
            "",
            f"   **Type**: {type_}",
            f"   **Default**: {default}",
        ]
        if unit:
            text.append(f"   **Unit**: {unit}")
        text.append("")
        text.append(f"   {description}")

        content = StringList(text)

        node = KeywordNode()
        self.state.nested_parse(content, self.content_offset, node)

        # Return both target and node to allow :ref: linking
        return [target_node, node]


def setup(app):
    app.add_node(
        KeywordNode,
        html=(visit_keyword_node_html, depart_keyword_node_html)
    )
    app.add_directive("keyword", KeywordDirective)
