from docutils import nodes
from docutils.parsers.rst import Directive, directives

class KeywordNode(nodes.General, nodes.Element): pass

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
        node = KeywordNode()
        text = f"**{keyword}**\n\n   *Type:* {type_}\n   {description}"
        self.state.nested_parse(
            self.state.input_lines.__class__([text]),
            self.content_offset,
            node
        )
        return [node]

def setup(app):
    app.add_node(KeywordNode, html=(visit_keyword_node_html, depart_keyword_node_html))
    app.add_directive("keyword", KeywordDirective)