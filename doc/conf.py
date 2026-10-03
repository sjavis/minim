# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Run doxygen -------------------------------------------------------------

import subprocess
subprocess.call("doxygen Doxyfile.in", shell=True)


# -- Project information -----------------------------------------------------

project = 'Minim'
copyright = '2021, Sam Avis'
author = 'Sam Avis'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['breathe', 'sphinx.ext.mathjax']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

highlight_language = 'c++'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Logo
html_logo = "logo/logo2_text.png"
html_theme_options = {'logo_only': True}


# -- Options for breathe -----------------------------------------------------

breathe_projects = {'minim': '_build/xml/'}
breathe_default_project = 'minim'
breathe_default_members = ('members',)


# -- Doxygen member groups as real sections ----------------------------------
# Breathe renders doxygen member-group (@name) titles as rubrics, which are
# not very prominent and do not appear in the page TOC. We extend breathe's
# class directive so that it emits the structure we want directly:
#
# - each class is wrapped in a real section (titled with the class name), and
# - each user-defined member group is rendered as a subsection of the class,
#
# giving prominent headings and correct nesting in the sidebar TOC.

def setup(app):
    from docutils import nodes
    from breathe.parser import DoxSectionKind
    from breathe.renderer.sphinxrenderer import SphinxRenderer

    # Render user-defined member groups (doxygen @name groups) as real
    # sections instead of rubrics, so they are prominent and appear in the
    # TOC nested under the class.
    def visit_sectiondef(self, node):
        node_list = []
        node_list.extend(self.render_optional(node.description))
        member_def = sorted(node.memberdef, key=lambda x: x.name) \
            if 'sort' in self.context.directive_args[2] else node.memberdef
        node_list.extend(self.render_iterable(member_def))
        if not node_list:
            return []
        if 'members-only' in self.context.directive_args[2]:
            return node_list
        text = self.section_titles[node.kind]
        if node.kind == DoxSectionKind.user_defined:
            text = node.header or 'Unnamed Group'
        idtext = text.replace(' ', '-').lower()
        if node.kind == DoxSectionKind.user_defined:
            # Prefix the id with the class name to keep anchors unique.
            # (get_qualification is empty when the class was found through
            # the doxygen index, so read the name from the context stack.)
            names = [
                n.value.compoundname for n in self.context.node_stack[1:]
                if hasattr(n, 'value') and hasattr(n.value, 'compoundname') and n.value.compoundname
            ]
            prefix = self.join_nested_name(names).replace('::', '-').lower()
            section = nodes.section(ids=[prefix + '-' + idtext])
            section += nodes.title('', text)
            # Wrap the group members in a nested definition list so they keep
            # the same styling as ungrouped members (the sphinx_rtd_theme
            # styles dls nested inside another dl differently to top-level
            # ones; members outside the class dl render blue instead of grey).
            # Any group description stays outside the wrapper.
            section.extend(node_list[:1])
            dl = nodes.definition_list()
            dl['classes'] = ['cpp', 'class']
            item = nodes.definition_list_item()
            definition = nodes.definition()
            definition.extend(node_list[1:])
            item += definition
            dl += item
            section += dl
            return [section]
        rubric = nodes.rubric(
            text=text,
            classes=['breathe-sectiondef-title'],
            ids=[idtext],
        )
        return [rubric] + node_list

    # Wrap the rendered class (a sphinx desc node) in a section titled with
    # the class name, so the group sections nest under the class in the TOC
    # and the class signature itself is the section heading.
    from breathe.directives.class_like import _DoxygenClassLikeDirective

    class SectionedClassDirective(_DoxygenClassLikeDirective):
        kind = 'class'

        def run(self):
            from sphinx import addnodes
            result = super().run()
            out = []
            for node_ in result:
                if not isinstance(node_, addnodes.desc):
                    out.append(node_)
                    continue
                sig = node_.next_node(addnodes.desc_signature)
                section = nodes.section()
                if sig is not None:
                    if sig['ids']:
                        section['ids'] = list(sig['ids'])
                        sig['ids'] = []
                    # Strip the namespace from the section title: the
                    # sidebar entry reads 'PhaseField' rather than
                    # 'minim::PhaseField'.
                    class_name = sig.get('_toc_name', '').split('::')[-1]
                    section += nodes.title('', class_name)
                node_['no-contents-entry'] = True
                section += node_
                # Move the user-defined group sections (rendered inside the
                # desc content, where sphinx's TOC collector skips them) out
                # into the class section, preserving document order.
                desc_content = node_.next_node(addnodes.desc_content)
                for group in list(desc_content.findall(nodes.section)):
                    if group.parent is not desc_content:
                        group.parent.remove(group)
                    else:
                        desc_content.remove(group)
                    section += group
                out.append(section)
            return out

    # Register our directive over breathe's (later registration wins).
    app.add_directive('doxygenclass', SectionedClassDirective, override=True)

    # Install the renderer override. breathe's @node_handler decorator fills
    # a dispatch registry (node_handlers) at class definition time, so the
    # method and the registry entry both need updating.
    from breathe.parser import Node_sectiondefType
    SphinxRenderer.visit_sectiondef = visit_sectiondef
    SphinxRenderer.node_handlers[Node_sectiondefType] = visit_sectiondef
