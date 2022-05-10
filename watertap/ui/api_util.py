"""
Utility functions for the ``api`` module.
"""
# standard library
import inspect
import json
import functools
from string import Template
from typing import Dict, Union, Optional, IO

# third-party
import fastjsonschema

#: Set this logger from the api module
util_logger = None


def log_meth(meth):
    @functools.wraps(meth)
    def wrapper(*args, **kwargs):
        name = _get_method_classname(meth)
        util_logger.debug(f"Begin {name}")
        try:
            result = meth(*args, **kwargs)
        except Exception as e:
            error_msg = f"Error in {name}"
            util_logger.exception(error_msg)
            raise e
        util_logger.debug(f"End {name}")
        return result

    return wrapper


def _get_method_classname(m):
    """Get class name for method, assuming method is bound and class has __dict__."""
    for k, v in inspect.getmembers(m):
        if k == "__qualname__":
            return v
    return "<unknown>"


# End logging


def open_file_or_stream(fos, attr, **kwargs) -> IO:
    """Open a file or use the existing stream. Avoids adding this logic to every function that wants to provide
       multiple ways of specifying a file.

    Args:
        fos: File or stream
        attr: Attribute to check on the ``fos`` object to see if it is a stream, e.g. "write" or "read"
        kwargs: Additional keywords passed to the ``open`` call. Ignored if the input is a stream.

    Returns:
        Opened stream object
    """
    if hasattr(fos, attr):
        output = fos
    else:
        output = open(fos, **kwargs)
    return output


# Utility function to automate documentation of something with a CONFIG that is not a ProcessBlock subclass


def config_docs(cls):
    """Class decorator to insert documentation for the accepted configuration options in the constructor's docstring.

    Returns:
        Decorator function
    """
    doc_lines = []
    tab_spc = " " * 4

    def format_list(item, depth):
        """Append to `doc_lines` during recursive walk of `item`."""
        if isinstance(item, tuple):
            indent = tab_spc * (depth + 3)
            bullet = ("-", "*")[depth == 1]
            name = item[0]
            desc = "(no description provided)" if len(item) == 1 else item[1]
            doc_lines.append(f"{indent}{bullet} `{name}`: {desc}")
        else:
            if depth > 0:
                doc_lines[-1] += " | items:"
            doc_lines.append("")
            for i in item:
                format_list(i, depth + 1)
            if depth > 0:
                doc_lines.append("")

    # Wrap in try/except so if something goes wrong no harm is done
    try:
        # Get Pyomo to generate 'documentation' that is parseable as a Python list
        s = cls.CONFIG.generate_documentation(
            indent_spacing=0,
            block_start="[",
            block_end="]",
            item_start="('%s',",
            item_body="'%s',",
            item_end="),",
        )
        # Parse the list then re-generate documentation as nested bulleted lists
        doc_list = eval(s)
        format_list(doc_list, 0)
    except Exception as err:
        util_logger.warning(f"Generating configuration docstring: {err}")
    # Add to the class constructor docstring. Assume that you can just append, i.e. the configuration
    # argument is the last entry in the "Args" section, which is the last section.
    if doc_lines:
        doc_str = "\n".join(doc_lines)
        cls.__init__.__doc__ = cls.__init__.__doc__.rstrip() + "\n" + doc_str
    # Return modified class
    return cls


# Schema validation
# -----------------


class SchemaException(Exception):
    """Exception for problems with schema definition."""

    pass


class JSONException(Exception):
    """Exception for problems encoding or decoding JSON."""

    pass


class Schema:
    """Lightweight wrapper for schema validation."""

    def __init__(self, definition: Union[Dict, str], **kwargs):
        """Constructor.

        Args:
            definition: Schema definition supported by underlying `fastjsonschema` parser. As of this writing,
                        this supports drafts 04, 06, and 07. This may be a string, which can be parsed as JSON,
                        or it can be a Python dict.
            kwargs: If present, ``Template.substutute()`` will be applied to the JSON-string version of the input
                    (which will be created if it starts as a Python dict), to substitute values in the
                    schema dynamically. This uses '$var' style substitution, with any unrecognized '$<thing>'
                    being left as-is.

        Raises:
            ValueError: Substitution of kwargs in schema fails
            JSONException: Parsing or formatting JSON fails (including converting to a string for applying the
                           substitution parameters).
            SchemaException: Compiling the schema fails
        """
        # Normalize input to Python dict, optionally after applying format params
        if kwargs:
            if hasattr(definition, "keys"):  # dict-like
                try:
                    definition_str = json.dumps(definition)
                except TypeError as err:
                    raise JSONException(
                        f"Converting input dict to JSON for formatting failed: {err}"
                    )
            else:
                definition_str = definition
            template = Template(definition_str)
            try:
                schema_str = template.safe_substitute(kwargs)
            except KeyError as err:
                raise ValueError(f"Substitution in schema definition failed: {err}")
            try:
                self._schema = json.loads(schema_str)
            except json.JSONDecodeError as err:
                raise JSONException(
                    f"Parsing of schema definition, after substitution, failed: {err}"
                )
        else:
            if hasattr(definition, "keys"):  # dict-like
                self._schema = definition
            else:
                try:
                    self._schema = json.loads(definition)
                except json.JSONDecodeError as err:
                    raise JSONException(
                        f"Parsing of schema definition after substitution of format params failed: {err}"
                    )
        # Compile input to a schema validation function
        try:
            self._validate = fastjsonschema.compile(self._schema)
        except fastjsonschema.JsonSchemaDefinitionException as err:
            raise SchemaException(f"Invalid schema definition: {err}")

    def validate(self, json_data: Union[Dict, str]) -> Optional[str]:
        """Validate input data against the schema given to the constructor.

        Args:
            json_data: Input in the form of a Python dict or JSON-format string.

        Returns:
            Validation error message, or None if it is valid (think of this as "no errors")

        Raises:
            JSONException: Input (string) data could not be decoded as JSON
        """
        if hasattr(json_data, "keys"):  # dict-like
            input_dict = json_data
        else:
            try:
                input_dict = json.loads(json_data)
            except json.JSONDecodeError as err:
                raise JSONException(f"Decoding input data failed: {err}")
        result = None
        try:
            self._validate(input_dict)
        except fastjsonschema.JsonSchemaValueException as err:
            result = err.message
        return result


# End schema validation
