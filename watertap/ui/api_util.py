"""
Utility functions for the ``api`` module.
"""
import inspect
import functools

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


