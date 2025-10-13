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
    option_spec = {'type': directives.unchanged}
    has_content = True

    def run(self):
        keyword = self.arguments[0]
        type_ = self.options.get('type', '')
        description = '\n'.join(self.content)

        # Build the reST content block that will be parsed into nodes
        text = [
            f"**{keyword}**",
            "",
            f"   *Type:* {type_}",
            f"   {description}",
        ]
        content = StringList(text)

        node = KeywordNode()
        self.state.nested_parse(content, self.content_offset, node)
        return [node]

def setup(app):
    app.add_node(
        KeywordNode,
        html=(visit_keyword_node_html, depart_keyword_node_html)
    )
    app.add_directive("keyword", KeywordDirective)
