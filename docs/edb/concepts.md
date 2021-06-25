# EDB Concepts

## Overview
The Electrolyte database stores metadata and data about chemical species, called here 
*components*, and *reactions*. It is accessed through a Python API
to return well-defined Python objects.

The data are stored in [MongoDB](https://mongodb.org), so they can be 
queried in a number of ways, and the system is extensible to new use-cases. The
native storage format for MongoDB is a [JSON](https://json.org) document, and the expected structure and
fields of the *component* and *reaction* data is defined by a [JSON Schema](https://json-schema.org). Validation
using those schemas is built into the API (though it can be disabled).

To interface with the [IDAES Core Modeling Framework](https://idaes-pse.readthedocs.io/en/stable/user_guide/concepts.html)
(IDAES-CMF, which underlies ProteusLib), add components and reactions to a "base" object
and fetch the result as a Python `dict`. This result that can be used to configure and 
build IDAES objects (ParameterBlocks, ReactionBlocks, etc.). The API also has methods to construct component and
reaction objects from IDAES configurations.

## Workflows
The EDB is intended to support some known workflows out of the box, with lower-level functions available when these
are not sufficient.

```{admonition} Coming soon
:class: note

This content is not yet finished.
```


