# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
#
# This module is a work in progress. Do not use it for real work right now.
#
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

"""
Electrolyte database.

Design::

     ┌────────────────────────────┐
     │           Validate         │
     │  ┌───────┐    ┌──────┐     │
     │  │JSON   │    │JSON  │     │
     │  │input  │    │schema│     │
     │  │data   │    │      │     │
     │  └───────┴┐ ┌─┴──────┘     │
     │           │ │              │
     │        ┌──▼─▼───┐          │
     │        │        │NO        │
     │        │ VALID? ├───►Error │
     │        └────┬───┘          │
     │             │              │
     │             │YES           │
     └─────────────┼──────────────┘
                   │                      ┌───────────────┐
     DB API........│.....                 │  reaction     │
        ...    ┌───▼───┐ ..    ;;;;;;;;   ├───────────────┤
       ..      │Load   ├──.───►;      ;   │  component    │
      ..       └───────┘  .    ; DB   ;───┼───────────────┤
     .          ▲         ..   ;      ;   │  base         │
    ..          │          .   ;;;;;;;;   ├───────────────┤
    .           │           ...           └───────────────┘
    .           │              ...        ...........................
    .   ┌───────▼┐    ┌─────────┐...    ...──────────┐   ┌───────────...
    .   │ Search ├───►│ Fetch   │  .   .. │ Component◄──o│ HasConfig │ .
    ..  └────────┘    └────┬────┘  .  ..  └──────▲───┘   └───o───────┘  .
     ...                   │      .  ..          │           │          .
        ...                │    ...  .           │           │         ..
           .....           │  ...   ..   ┌───────x───┐       │         .
                 . ... .. .│...     .    │ Reaction  │◄──────┘        ..
                           │        .    └─────▲─────┘                .
                           │        .          │                    ..
                           │        .  **┌─────x───┐               ..
                           └─────────.──►│ Result  │     Data model
                                     ... └─────────┘            ..
                                       .....                ....
                                            ............ ...
"""
