"""Render author-year citations with round brackets: Sims (1972), (Sims, 1972).

sphinxcontrib-bibtex's built-in ``author_year`` style uses square brackets.
The book's prose has always written citations as "Sims (1972)", so we register
a variant that differs from ``author_year`` only in the bracket characters.
"""
from dataclasses import dataclass, field

import sphinxcontrib.bibtex.plugin
from sphinxcontrib.bibtex.style.referencing import BracketStyle
from sphinxcontrib.bibtex.style.referencing.author_year import (
    AuthorYearReferenceStyle,
)


def parens() -> BracketStyle:
    return BracketStyle(left="(", right=")")


@dataclass
class AuthorYearParenReferenceStyle(AuthorYearReferenceStyle):
    bracket_textual: BracketStyle = field(default_factory=parens)
    bracket_parenthetical: BracketStyle = field(default_factory=parens)
    bracket_author: BracketStyle = field(default_factory=parens)
    bracket_label: BracketStyle = field(default_factory=parens)
    bracket_year: BracketStyle = field(default_factory=parens)


# Register at import time: the style is resolved when the builder initialises.
sphinxcontrib.bibtex.plugin.register_plugin(
    "sphinxcontrib.bibtex.style.referencing",
    "author_year_paren",
    AuthorYearParenReferenceStyle,
)


def setup(app):
    return {"parallel_read_safe": True, "parallel_write_safe": True}
