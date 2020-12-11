"""
NAWI Electrolyte Database
"""
import argparse
from collections import namedtuple
import json
import logging
from pathlib import Path
import random
import sqlite3
import time
from typing import Union, Sequence, List, Dict

__author__ = "Dan Gunter"

# Type aliases
Filename = Union[Path, str]

_log = logging.getLogger(__name__)


class DuplicateReactionError(Exception):
    def __init__(self, existing_id, new_id):
        msg = f"Reaction id={existing_id} has the same components as reaction={new_id}"
        super().__init__(msg)


class VersionMismatch(Exception):
    def __init__(self, code, db):
        msg = f"Database version {db} conflicts with code version {code}"
        super().__init__(msg)


class Component:
    """A reaction component (chemical species).
    """

    def __init__(
        self,
        formula,
        id_=None,
        num=1,
        type_=None,
        casrn: str = None,
        molecular_weight: float = None,
        charge: int = None,
        dissociation: str = None,
    ):
        self.formula = formula
        self.id_ = id_
        self.type = type_
        self.casrn = casrn
        self.molecular_weight = molecular_weight
        self.charge = charge
        self.dissociation = dissociation
        self.num = num

    def __str__(self):
        if self.num == 1:
            return f"{self.formula}"
        else:
            return f"({self.formula}){self.num}"

class Reaction:
    """Representation of a reaction and associated properties.
    """

    ACID_BASE_VALUES = ("acidic", "basic", "any")

    def __init__(
        self,
        electrolyte_db: "ElectrolyteDB",
        id_=None,
        description="",
        form: str = "general",
        components = None,
        conditions: str = "",
        acid_base: str = ACID_BASE_VALUES[2],
        pka298k: float = 0.0,
    ):
        assert (
            acid_base in self.ACID_BASE_VALUES
        ), "Invalid value for `acid_base` argument"
        self._db = electrolyte_db
        self.id_ = id_
        self.properties = {
            "form": form,
            "description": description,
            "conditions": conditions,
            "acid_base": acid_base,
            "pka298k": pka298k,
        }
        # Accept a list of ids or components (or None) for
        # the 'component' parameter
        if components is None or len(components) == 0:
            self._comp, self._comp_ids = None, None
        elif isinstance(list(components)[0], int):
            self._comp, self._comp_ids = None, frozenset(components)
        else:
            self._comp, self._comp_ids = frozenset(components), None
        self._details = None  # lazily initialized with reaction details
        self._lhs, self._rhs = None, None

    @property
    def component_ids(self):
        if self._comp_ids is None:
            _ = self.components  # this sets self._comp_ids
        return self._comp_ids

    @property
    def components(self):
        if self._comp is None:
            self._fetch()
            component_list, formulae = [], set()
            for side_name in "lhs", "rhs":
                for comp_detail in self._details[side_name]:
                    if comp_detail["formula"] in formulae:
                        continue  # skip duplicates
                    component_list.append(Component(**comp_detail))
                    formulae.add(comp_detail["formula"])
            self._comp = frozenset(component_list)
            self._comp_ids = frozenset([x.id_ for x in component_list])
        return self._comp

    @property
    def formula(self) -> (Dict, Dict):
        return self.lhs, self.rhs

    @property
    def lhs(self) -> Dict:
        if self._lhs is None:
            self._fetch()
            self._lhs = self._extract_formula(self._details["lhs"])
        return self._lhs

    @property
    def rhs(self) -> Dict:
        if self._rhs is None:
            self._fetch()
            self._rhs = self._extract_formula(self._details["rhs"])
        return self._rhs

    @staticmethod
    def _extract_formula(details_list):
        return {item["formula"]: item["num"] for item in details_list}

    def _fetch(self):
        if self._details is None:
            self._details = self._db.fetch_reaction_details(self.id_, self.properties["form"])

    def as_html(self, properties=False):
        return ""  # TODO: Return an HTML version of the formula

    def as_latex(self, properties=False):
        return ""  # TODO: Return a LaTeX version of the formula

    def as_text(self, properties=False):
        formula_text = " <-> ".join((" + ".join(self.lhs), " + ".join(self.rhs)))
        lines = [
            ("ID", self.id_),
            ("Description", self.properties["description"]),
            ("Formula", formula_text),
            ("Form", self.properties["form"]),
            ("Conditions", self.properties["conditions"]),
            ("Acid/Base", self.properties["acid_base"]),
            ("pKa at 298K", self.properties["pka298k"])
        ]
        return self._twocol(lines, col_sep=": ")

    @staticmethod
    def _twocol(lines, col_sep="", row_sep="\n"):
        w = max((len(s[0]) for s in lines))
        rows = [("{0:<%d}{sep}{1}" % w).format(sep=col_sep, *keyval) for keyval in lines]
        return row_sep.join(rows)

    def __str__(self):
        return self.as_text()


class VersionNumber(namedtuple("VersionNumber", ["major", "minor"])):
    """Simple version number.
    """

    def __str__(self):
        return f"{self.major}.{self.minor}.0"  # add dummy PATCH version


class ElectrolyteDB:
    """Methods for accessing the electrolyte database.
    """

    VERSION = VersionNumber(0, 1)

    DB_LIST_SEP = ";"  # use this to scrunch lists together as strings

    def __init__(self, file: Filename):
        """Connect to database.

        The database will be created and initialized if it does not exist.

        Args:
            file: SQLite3 database file

        Raises:
             VersionMismatch, if major versions do not match between this code and an existing DB.
        """
        self._path = Path(file)
        db_exists = self._path.exists()
        self._db = sqlite3.connect(self._path)
        self._comp_ids, self._comp_fml = {}, {}
        if db_exists:
            self._check_version()
        else:
            self._init_tables()
        self._fetch_component_ids()
        self._comp_react = {}  # {<side>: {<form>: {component-id-set: [reaction ids.. ]} } }
        self._needs_update = True

    def get_component(self, component_formula) -> Union[Component, None]:
        """Lookup component in database.

        Args:
            component_formula: Component formula text, e.g., "H2O" or "KHCO3"

        Returns:
            A component object, or None if no match is found in database
        """
        self._update()
        cursor = self._db.cursor()
        cursor.execute(
            "SELECT id, type, dissociation, charge, casrn, mwt FROM component WHERE fml=:fml", {"fml": component_formula}
        )
        row = cursor.fetchone()
        if row is None:
            return None
        id_, type_, diss, charge, casrn, mwt = row
        diss_list = [] if diss is None else diss.split(self.DB_LIST_SEP)
        return Component(component_formula, id_=id_, type_=type_, molecular_weight=mwt, casrn=casrn,
                         charge=charge, dissociation=diss_list)

    def get_reactions(
        self, components: Sequence[Union[Component, str]], form_filter=None, match_all: bool = False
    ) -> List[Reaction]:
        """Get reactions based on a component list.

        Args:
            components: Find reactions for which all components are in this list
            form_filter: Iterable of allowed forms of reaction, or None/empty for any
            match_all: If true, components must include all components from both sides of the reaction.
                       Otherwise, components must only match all components on either side (the default).

        Returns:
            The reactions matching the components.
        """
        self._update()
        # normalize value for 'form_filter' to either None or a set()
        if form_filter is not None:
            if len(form_filter) == 0:
                form_filter = None
            else:
                form_filter = set(form_filter)
        # In input components, replace anything that doesn't have an `id_` attribute with Component object
        # looked up from the str() of its value. Usually this will be just a string.
        component_set = set()
        for c in components:
            if isinstance(c, int):
                id_ = c
            elif not hasattr(c, "id_"):
                id_ = self.get_component(str(c)).id_
            else:
                id_ = c.id_
            component_set.add(id_)
        component_key = frozenset(component_set)
        # Iterate over all forms
        found_reactions_all_forms = {}
        if match_all:
            # only the 'both' will match
            for form, reaction_list in self._find_reactions(component_key, self._comp_react['both']):
                found_reactions_all_forms[form] = reaction_list
        else:
            # either LHS or RHS will match
            for side in 'lhs', 'rhs':
                # print(f"@@ {side} match")
                for form, reaction_list in self._find_reactions(component_key, self._comp_react[side]):
                    found_reactions_all_forms[form] = reaction_list
        reactions = []
        for form, reaction_ids in found_reactions_all_forms.items():
            if form_filter is not None and form not in form_filter:
                # wrong form; skip
                continue
            cursor = self._db.cursor()
            cursor.execute(f"SELECT r.description, r.conditions, r.acid_base, r.pka298k FROM reaction r "
                           f"WHERE r.id IN ({','.join([str(rid) for rid, _ in reaction_ids])})")
            prop_kwds = []  # array of dict; keys matching keywords in Reaction constructor
            for row in cursor.fetchall():
                prop_kwds.append({"description": row[0], "conditions": row[1], "acid_base": row[2], "pka298k": row[3]})
            for prop_kwd, (r_id, c_id) in zip(prop_kwds, reaction_ids):
                r = Reaction(self, id_=r_id, form=form, components=c_id, **prop_kwd)
                reactions.append(r)
        # Done
        return reactions

    @staticmethod
    def _find_reactions(component_key, cr):
        """Find reactions where components are a subset of the component_key set.
        """
        for form in cr:
            found_reactions = []
            for component_ids, reaction_id in cr[form].items():
                if component_ids.issubset(component_key):
                    # print(f"@@ {component_ids} is-subset {component_key}: r = {reaction_id}")
                    try:
                        # test 'reaction_id' for iterable (duck typing)
                        iter(reaction_id)
                    except TypeError:
                        # not iterable
                        found_reactions.append((reaction_id, component_ids))
                    else:
                        # iterable
                        for rid in reaction_id:
                            found_reactions.append((rid, component_ids))
            # add if found
            if found_reactions:
                yield form, found_reactions

    def fetch_reaction_details(self, reaction_id, form):
        """Get details of components in the reaction.
        """
        stmt = (
            f"SELECT c.fml, c.casrn, c.mwt, c.type, c.charge, c.dissociation, "
            f"cs.comp_num, cs.side, c.id "
            f"FROM reaction r JOIN reaction_form rf ON (r.id = rf.react_id) "
            f"JOIN form f ON (rf.form_id = f.id and f.form = '{form}') "
            f"JOIN component_side cs ON (f.lhs_id = cs.id OR f.rhs_id = cs.id) "
            f"JOIN component c ON (cs.comp_id = c.id) "
            f"WHERE r.id = {reaction_id}"
        )
        c = self._db.cursor()
        c.execute(stmt)
        rows = c.fetchall()
        if rows is None:
            result = None
        else:
            result = {"lhs": [], "rhs": []}
            for row in rows:
                # These keywords should match attribute names in Component class,
                # since this fact is used to construct a Component from this dict
                # TODO: Just make this a Component?
                component_details = {
                    "formula": row[0],
                    "casrn": row[1],
                    "molecular_weight": row[2],
                    "type_": row[3],
                    "charge": row[4],
                    "dissociation": row[5],
                    "num": row[6],
                    "id_": row[7],
                }
                side = row[7]
                assert side in (0, 1), "Side must be 0 (zero) or 1"
                (result["lhs"], result["rhs"])[side].append(component_details)
        return result

    def is_empty(self):
        self._update()
        # add up length for each 'side' sub-dict
        return sum(map(len, self._comp_react.values())) == 0

    def _init_tables(self):
        """Build the schema for a new database.
        """
        script = f"""
            CREATE TABLE version (major INT, minor INT);
            CREATE TABLE component (id INTEGER PRIMARY KEY AUTOINCREMENT, fml TEXT, casrn TEXT, mwt FLOAT, type TEXT, 
                charge INT, dissociation TEXT);
            CREATE TABLE component_side (id INTEGER, comp_id INT, comp_num INT, side INT, FOREIGN KEY(comp_id)
                REFERENCES component(id));
            CREATE TABLE form (id INTEGER PRIMARY KEY AUTOINCREMENT, lhs_id INT, rhs_id INT, form TEXT, 
                FOREIGN KEY(lhs_id) REFERENCES component_side(id),FOREIGN KEY(rhs_id) REFERENCES component_side(id));
            CREATE TABLE reaction (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, description TEXT, conditions TEXT,
                acid_base TEXT, pka298k FLOAT);
            CREATE TABLE reaction_form (react_id INT, form_id INT, FOREIGN KEY(react_id) REFERENCES reaction(id), 
                FOREIGN KEY(form_id) REFERENCES form(id));
            CREATE INDEX cs1 ON component_side (id);
            CREATE INDEX form_lhs_idx ON form (lhs_id);
            CREATE INDEX form_rhs_idx ON form (rhs_id);
            CREATE INDEX reaction_form_idx ON reaction_form (react_id, form_id);
            INSERT INTO version VALUES ({self.VERSION.major}, {self.VERSION.minor});
        """
        _log.debug("Creating schema in new database")
        self._db.executescript(script)
        _log.debug("Database schema created")

    def _check_version(self):
        c = self._db.cursor()
        c.execute("SELECT major, minor FROM version")
        row = c.fetchone()
        db_version = VersionNumber(major=row[0], minor=row[1])
        if self.VERSION.major != db_version.major:
            raise VersionMismatch(self.VERSION, db_version)
        if self.VERSION.minor < db_version.minor:
            _log.warning(
                f"Database version {db_version} is newer than code version {self.VERSION}"
            )
        elif self.VERSION.minor > db_version.minor:
            _log.warning(
                f"Database version {db_version} is older than code version {self.VERSION}"
            )

    def load_from_file(self, file: Filename):
        path = Path(file)
        data = json.load(path.open())
        self.load_from_json(data)

    def load_from_json(self, data: dict):
        if "components" in data:
            self._load_components(data["components"])
        if "reactions" in data:
            self._load_reactions(data["reactions"], data["components"])
        self._needs_update = True

    def _load_components(self, comp):
        values = []
        for c in comp:
            row = [
                None,
                c["name"],
                c["CASRN"],
                c["molecular_weight"],
                c["type"],
                c["charge"],
            ]
            if c["dissociation"] is None:
                row.append(None)
            else:
                row.append(self.DB_LIST_SEP.join(c["dissociation"]))
            values.append(row)
        self._db.executemany("INSERT INTO component VALUES (?,?,?,?,?,?,?)", values)
        self._db.commit()

    def _load_reactions(self, react, comp):
        self._fetch_component_ids()  # populates self._comp_ids
        cursor = self._db.cursor()
        # get starting reaction_side id
        cursor.execute("SELECT MAX(id) from component_side")
        max_id = cursor.fetchone()[0]
        rs_id = 1 if max_id is None else max_id + 1
        # add all the reactions to both the 'component_side' and 'reaction' tables
        for r in react:
            rform_ids = []  # remember ids of rows added to 'form' table
            for form in r["forms"]:
                side_ids = [0, 0]
                for side in r["forms"][form]:
                    formulae = r["forms"][form][side]
                    side_int = 0 if side == "LHS" else 1
                    for fml, num in formulae.items():
                        comp_id = self._comp_ids[fml]
                        values = [rs_id, comp_id, num, side_int]
                        cursor.execute(
                            "INSERT INTO component_side VALUES (?,?,?,?)", values
                        )
                    side_ids[side_int] = rs_id
                    rs_id += 1  # new ID for other side (or next reaction)
                values = [None] + side_ids + [form]
                cursor.execute("INSERT INTO form VALUES (?,?,?,?)", values)
                rform_ids.append(cursor.lastrowid)
            # add a new reaction, and get its auto-generated identifier
            cursor.execute(
                "INSERT INTO reaction VALUES (?, ?, ?, ?, ?, ?)",
                [None, r["ID"], r["description"], r["conditions"], r["acid_base"], r["pKa_298K"]],
            )
            rct_id = cursor.lastrowid
            # add entries to join table for reaction and its forms
            for id_ in rform_ids:
                cursor.execute("INSERT INTO reaction_form VALUES (?,?)", [rct_id, id_])
        cursor.close()
        self._db.commit()

    def _fetch_component_ids(self):
        self._comp_ids, self._comp_fml = {}, {}
        for row in self._db.execute("SELECT id, fml FROM component"):
            self._comp_ids[row[1]] = row[0]
            self._comp_fml[row[0]] = row[1]

    def _update(self):
        if not self._needs_update:
            return
        self._comp_react = {'lhs': {}, 'rhs': {}, 'both': {}}
        for side in self._comp_react:
            # Query for LHS only, RHS only, or both is identical except this part
            if side == "lhs":
                side_query = "(f.lhs_id == s.id)"
                unique_reactions = False  # may be multiple reactions for same LHS
            elif side == "rhs":
                side_query = "(f.rhs_id == s.id)"
                unique_reactions = False  # may be multiple reactions for same RHS
            else:
                side_query = "(f.lhs_id == s.id or f.rhs_id == s.id)"
                unique_reactions = True  # may NOT be multiple reactions for same LHS + RHS
            cursor = self._db.cursor()
            # Join the 'form' and 'component_side' tables to get components for each reaction-form,
            # and then use the 'reaction_form' join table to get the (parent) reaction identifier
            _log.debug("Getting reaction+form+component with a database join..")
            t0 = time.time()
            cursor.execute(
                f"SELECT rf.react_id, f.form, s.comp_id "
                f"FROM form AS f "
                f"JOIN component_side AS s ON {side_query} "
                f"JOIN reaction_form as rf ON (rf.form_id == f.id) "
                f"ORDER BY rf.react_id"
            )
            t1 = time.time()
            _log.debug(f"Join for side={side} took {t1 - t0:.1f}s")
            _log.debug("Processing reactions..")
            t0 = time.time()
            prev_react_id, forms, comp_ids = None, None, None
            # Accumulate components for each reaction, and each 'form' of a reaction
            for react_id, react_form, comp_id in cursor.fetchall():
                # if this form has not been encountered, make a new dict for it
                # (this should be fast: only a few keys to check)
                if react_form not in self._comp_react[side]:
                    self._comp_react[side][react_form] = {}
                # When we encounter a new reaction, record the previous one
                if react_id != prev_react_id:
                    # If there is a previous reaction to record..
                    if prev_react_id is not None:
                        self._update_add_reaction_components(prev_react_id, comp_ids, self._comp_react[side],
                                                             unique=unique_reactions)
                    prev_react_id, comp_ids = react_id, {}
                # Add components from this reaction
                if react_form not in comp_ids:
                    comp_ids[react_form] = {comp_id}
                else:
                    comp_ids[react_form].add(comp_id)
            t1 = time.time()
            _log.debug(f"Processing reactions for side={side} took {t1 - t0:.1f}s")
        self._needs_update = False

    @staticmethod
    def _update_add_reaction_components(react_id, comp_ids, comp_react, unique=False):
        """For _update(), add all the accumulated reaction components for each 'form'
        """
        for reaction_form, component_ids in comp_ids.items():
            frz_component_ids = frozenset(component_ids)
            comp_react_form = comp_react[reaction_form]
            if frz_component_ids in comp_react_form:
                if unique:
                    other_react_id = comp_react_form[frz_component_ids]
                    raise DuplicateReactionError(react_id, other_react_id)
            elif not unique:
                comp_react_form[frz_component_ids] = []
            if unique:
                comp_react_form[frz_component_ids] = react_id
            else:
                comp_react_form[frz_component_ids].append(react_id)


def __test(dbfile, inputfile, stress=False, delete=False):
    if delete:
        _log.debug(f"Deleting existing DB file '{dbfile}'")
        p = Path(dbfile)
        p.unlink(missing_ok=True)
    edb = ElectrolyteDB(dbfile)
    if stress:
        # stress-test
        component_list = []
        letters = "ABCDEFGHIJK"

        def make_formula(i):
            return "".join([letters[int(c)] for c in str(i)])

        num_components = 10000
        for i in range(num_components):
            c = {
                "CASRN": "7732-18-5",
                "name": make_formula(i),
                "alternate_names": [],
                "molecular_weight": 18.01528,
                "type": random.choice(("molecular", "anion", "cation")),
                "charge": random.choice([None, 1, -1]),
                "dissociation": None,
            }
            component_list.append(c)
        reaction_list = []
        reaction_components = set()
        num_reactions = 10000
        for rnum in range(num_reactions):
            while True:  # until unique set of reaction components chosen
                chosen_components = set()
                comp = []
                for i in range(5):
                    while True:
                        j = random.randint(0, num_components - 1)
                        if j not in chosen_components:
                            break
                    chosen_components.add(j)
                    comp.append(make_formula(j))
                comp = tuple(sorted(comp))
                if comp not in reaction_components:
                    break
            reaction_components.add(comp)
            r = {
                "ID": f"L{rnum:04}",
                "description": "Self ionization of nothing",
                "type": "liquid phase",
                "forms": {
                    "hydronium": {"LHS": {comp[0]: 2}, "RHS": {comp[1]: 1, comp[2]: 1}},
                    "proton": {"LHS": {comp[0]: 1}, "RHS": {comp[3]: 1, comp[4]: 1}},
                },
                "conditions": "equilibrium reaction upon dissolution of CO2 in water",
                "acid_base": "any",
                "pKa_298K": None,
            }
            reaction_list.append(r)
        t0 = time.time()
        edb.load_from_json({"components": component_list, "reactions": reaction_list})
        t1 = time.time()
        delta_t = t1 - t0
        print(
            f"Time to load {num_components} components and {num_reactions} reactions: {delta_t:.1f}s"
        )
        # query
        # 0. update
        print("Running update")
        t0 = time.time()
        edb.get_component("DEADBEEF")
        t1 = time.time()
        delta_t = t1 - t0
        print(
            f"Time to update with {num_components} components and {num_reactions} reactions: {delta_t:.1f}s"
        )
        # 1. get_component
        t0 = time.time()
        for i in range(100):
            cnum = random.randint(0, num_components - 1)
            cfml = make_formula(cnum)
            edb.get_component(cfml)
        delta_t = time.time() - t0
        print(f"Time to fetch 100 components: {delta_t:.1f}s")
        # 2. get_reaction IDs
        reps = 100
        print(f"Getting reaction IDs {reps} times..", end="")
        total_time = 0
        all_reaction_ids = []
        for _ in range(reps):
            comp_fml, form = set(), None
            for i in range(3):
                r = random.choice(reaction_list)
                # print(f"Choose reaction: {r}")
                if form is None:
                    form = random.choice(list(r["forms"].keys()))
                rf = r["forms"][form]
                for side in "LHS", "RHS":
                    for k in rf[side]:
                        comp_fml.add(k)
            comp_fml_str = ",".join(comp_fml)
            comp = [edb.get_component(fml) for fml in comp_fml]
            t0 = time.time()
            reactions = edb.get_reactions(comp)
            total_time += time.time() - t0
            # print(f"Found {n} reactions")
            print(".", end="", flush=True)
        print(f"\nTime to fetch reaction IDs {reps} times: {total_time:.1f}s")
        # 3. Get first 50 and last 50 reactions
        n = 50
        print(f"Get reaction details {n * 2} times", end="")
        t0 = time.time()
        for i in range(1, n + 1):
            react = Reaction(edb, id_=i)
            fml = react.formula
            print(".", end="", flush=True)
        t1 = time.time()
        for i in range(num_reactions, num_reactions - n, -1):
            react = Reaction(edb, id_=i)
            fml = react.formula
            print(".", end="", flush=True)
        t2 = time.time()
        # print(f"Final formula: {react.formula}")
        print(
            f"\nTime to get first {n} reactions = {t1 - t0:.4f}s and last {n} reactions = {t2 - t1:.4f}s"
        )
        print(f"Total time to get {n * 2} reactions = {t2 - t0:.4f}s")
    else:
        # non-stress test
        if inputfile:
            edb.load_from_file(inputfile)
        elif delete:
            print("ERROR: deleted database, but no input file provided")
            return
        if edb.is_empty():
            print("ERROR: Database is empty!")
            return
        components_i_like = ("OH-", "H2O", "H3O+", "CO2", "HCO3-")
        components = []
        for c in components_i_like:
            component = edb.get_component(c)
            print(f"Component {c} type={component.type}")
            components.append(component)
        rr = edb.get_reactions(components)
        print(f"Reactions with {components_i_like}: {', '.join(map(str, rr))}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("dbfile")
    ap.add_argument("inputfile", nargs="?", default=None)
    ap.add_argument("--stress", action="store_true")
    ap.add_argument("--delete", "-d", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0, dest="vb")
    args = ap.parse_args()
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    _log.addHandler(_h)
    if args.vb > 1:
        _log.setLevel(logging.DEBUG)
    elif args.vb > 0:
        _log.setLevel(logging.INFO)
    else:
        _log.setLevel(logging.WARNING)
    __test(args.dbfile, args.inputfile, stress=args.stress, delete=args.delete)
