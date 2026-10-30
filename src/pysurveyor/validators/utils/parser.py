# Coding: UTF-8

# Licensed under the GNU General Public License v3.0.
# Full text of the license can be found in the LICENSE and COPYING files.

from __future__ import annotations

import ast
from collections.abc import Callable
import operator
import re
from typing import Any

from pysurveyor.backend.backend import NA
from pysurveyor.exceptions.exceptions import ParseError
from pysurveyor.validators.abc.validator import Validator
from pysurveyor.validators.between import Between
from pysurveyor.validators.comparison_validators import (
    Equal,
    Greater,
    GreaterOrEqual,
    Less,
    LessOrEqual,
)
from pysurveyor.validators.conversion_validators import (
    IsBoolean,
    IsNumeric,
    IsText,
)
from pysurveyor.validators.utils.constants import COMMENT

# ======================================================================================
# Available validators
# ======================================================================================
_VALIDATORS = {
    cls.__name__: cls
    for cls in (
        IsNumeric,
        IsBoolean,
        IsText,
        Between,
        Equal,
        Greater,
        GreaterOrEqual,
        Less,
        LessOrEqual,
    )
}


# ======================================================================================
# Supported operators
# ======================================================================================
_BINARY_OPERATORS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.FloorDiv: operator.floordiv,
    ast.Mod: operator.mod,
    ast.Pow: operator.pow,
    ast.BitAnd: operator.and_,
}

_UNARY_OPERATORS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
    ast.Not: operator.not_,
}

_COMPARISON_OPERATORS = {
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
    ast.In: lambda left, right: left in right,
    ast.NotIn: lambda left, right: left not in right,
    ast.Is: operator.is_,
    ast.IsNot: operator.is_not,
}

_LAMBDA_NODE_TYPES = (
    ast.Expression,
    ast.Lambda,
    ast.arguments,
    ast.arg,
    ast.Constant,
    ast.Name,
    ast.Attribute,
    ast.List,
    ast.Tuple,
    ast.Set,
    ast.Dict,
    ast.BinOp,
    ast.UnaryOp,
    ast.BoolOp,
    ast.Compare,
    ast.Load,
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.FloorDiv,
    ast.Mod,
    ast.Pow,
    ast.BitAnd,
    ast.UAdd,
    ast.USub,
    ast.Not,
    ast.And,
    ast.Or,
    ast.Eq,
    ast.NotEq,
    ast.Lt,
    ast.LtE,
    ast.Gt,
    ast.GtE,
    ast.In,
    ast.NotIn,
    ast.Is,
    ast.IsNot,
)


# ======================================================================================
# Restricted expression evaluator
# ======================================================================================
class _Evaluator:

    def evaluate(
        self,
        node: ast.AST,
        scope: dict[str, Any] | None = None,
    ) -> Any:

        scope = scope or {}

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Constant):
            value = node.value
            if value is Ellipsis:
                raise ParseError("Ellipsis is not supported")
            return value

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Name):
            if node.id in scope:
                return scope[node.id]
            if node.id == "None":
                return None
            if node.id == "True":
                return True
            if node.id == "False":
                return False
            if node.id == "COMMENT":
                return COMMENT
            raise ParseError(f"Unknown name: {node.id!r}")

        if isinstance(node, ast.Attribute):
            # backend missing-value are recognized -------------------------------------
            if isinstance(node.value, ast.Name) and node.value.id in ["pd", "pandas"]:
                if node.attr == "NA":
                    return NA

            # regex flags are recognized -----------------------------------------------
            if isinstance(node.value, ast.Name) and node.value.id == "re":
                regex_flags = {
                    "IGNORECASE": re.IGNORECASE,
                    "MULTILINE": re.MULTILINE,
                    "DOTALL": re.DOTALL,
                    "VERBOSE": re.VERBOSE,
                    "ASCII": re.ASCII,
                    "UNICODE": re.UNICODE,
                }
                if node.attr in regex_flags:
                    return regex_flags[node.attr]

            raise ParseError("Attribute access is not allowed")

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.List):
            return [self.evaluate(item, scope) for item in node.elts]

        if isinstance(node, ast.Tuple):
            return tuple(self.evaluate(item, scope) for item in node.elts)

        if isinstance(node, ast.Set):
            return {self.evaluate(item, scope) for item in node.elts}

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Dict):
            return {
                self.evaluate(key, scope): self.evaluate(value, scope)
                for key, value in zip(node.keys, node.values)
            }

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Lambda):
            return self._lambda(node, scope)

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Call):
            return self._call(node, scope)

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.BinOp):
            operation = _BINARY_OPERATORS.get(type(node.op))

            if operation is None:
                raise ParseError("Unsupported binary operator")

            left = self.evaluate(node.left, scope)
            right = self.evaluate(node.right, scope)

            try:
                return operation(left, right)
            except Exception as e:
                raise ParseError(f"Error evaluating binary operation: {e}") from e

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.UnaryOp):
            operation = _UNARY_OPERATORS.get(type(node.op))

            if operation is None:
                raise ParseError("Unsupported unary operator")

            try:
                return operation(self.evaluate(node.operand, scope))
            except Exception as e:
                raise ParseError(f"Error evaluating unary operation: {e}") from e

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.BoolOp):

            if isinstance(node.op, ast.And):
                result: Any = True

                for value in node.values:
                    result = self.evaluate(value, scope)

                    if not result:
                        return result

                return result

            if isinstance(node.op, ast.Or):
                result = False

                for value in node.values:
                    result = self.evaluate(value, scope)

                    if result:
                        return result

                return result

            raise ParseError("Unsupported boolean operator")

        # ------------------------------------------------------------------------------
        if isinstance(node, ast.Compare):
            left = self.evaluate(node.left, scope)

            for op_node, comparator in zip(node.ops, node.comparators):

                right = self.evaluate(comparator, scope)
                operation = _COMPARISON_OPERATORS.get(type(op_node))

                if operation is None:
                    raise ParseError("Unsupported comparison operator")

                if not operation(left, right):
                    return False

                left = right

            return True

        raise ParseError(f"Unsupported expression: {type(node).__name__}")

    # ==================================================================================
    def _lambda(
        self,
        node: ast.Lambda,
        parent_scope: dict[str, Any],
    ) -> Callable[..., Any]:

        args = node.args

        # ------------------------------------------------------------------------------
        if (
            args.posonlyargs
            or args.vararg is not None
            or args.kwarg is not None
            or args.kwonlyargs
            or args.defaults
            or args.kw_defaults
        ):
            raise ParseError(
                "Lambdas must use simple positional arguments without defaults"
            )

        names = [arg.arg for arg in args.args]

        if len(names) != len(set(names)):
            raise ParseError("Duplicate lambda argument")

        # ------------------------------------------------------------------------------
        for child in ast.walk(node.body):

            if not isinstance(child, _LAMBDA_NODE_TYPES):
                raise ParseError(
                    f"Unsupported lambda operation: {type(child).__name__}"
                )

            if isinstance(child, ast.Name) and child.id not in names:
                raise ParseError(f"Unknown name in lambda: {child.id!r}")

        # ------------------------------------------------------------------------------
        def predicate(*values: Any) -> Any:
            if len(values) != len(names):

                raise TypeError(
                    f"Lambda expects {len(names)} argument(s), "
                    f"received {len(values)}"
                )

            local_scope = dict(parent_scope)
            local_scope.update(zip(names, values))

            return self.evaluate(node.body, local_scope)

        # ------------------------------------------------------------------------------
        return predicate

    # ==================================================================================
    def _call(self, node: ast.Call, scope: dict[str, Any]) -> Any:

        if any(keyword.arg is None for keyword in node.keywords):
            raise ParseError("Argument unpacking is not supported")

        # ------------------------------------------------------------------------------
        if isinstance(node.func, ast.Name):
            name = node.func.id
            if name not in _VALIDATORS:
                raise ParseError(f"Unknown or forbidden function: {name!r}")
            function: Callable[..., Any] = _VALIDATORS[name]

        elif (
            isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "re"
            and node.func.attr == "compile"
        ):
            function = re.compile

        else:
            raise ParseError("Function call is not allowed")

        # ------------------------------------------------------------------------------
        args = [self.evaluate(arg, scope) for arg in node.args]
        kwargs = {
            keyword.arg: self.evaluate(keyword.value, scope)
            for keyword in node.keywords
        }

        # ------------------------------------------------------------------------------
        try:
            result = function(*args, **kwargs)
        except Exception as exc:
            raise ParseError(
                f"Could not construct {getattr(function, '__name__', function)!s}: "
                f"{exc}"
            ) from exc

        # ------------------------------------------------------------------------------
        return result


# ======================================================================================
# Public API
# ======================================================================================
def parse_validator(expression: str) -> Validator:
    """Parse a string and return the corresponding Validator object."""

    expression = expression.strip()

    # ----------------------------------------------------------------------------------
    if not isinstance(expression, str):
        raise TypeError("expression must be a string")
    if not expression:
        raise ParseError("Expression cannot be empty")

    # ----------------------------------------------------------------------------------
    try:
        tree = ast.parse(expression, mode="eval")
    except (SyntaxError, ValueError) as exc:
        raise ParseError(f"Invalid expression: {exc}") from exc

    # ----------------------------------------------------------------------------------
    evaluator = _Evaluator()
    result = evaluator.evaluate(tree.body)

    # ----------------------------------------------------------------------------------
    if not isinstance(result, Validator):
        raise ParseError("Expression must produce a Validator instance")

    # ----------------------------------------------------------------------------------
    return result
