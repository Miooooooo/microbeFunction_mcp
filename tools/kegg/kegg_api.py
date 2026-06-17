"""KEGG REST API wrapper."""

from __future__ import annotations

import time

import requests


class KeggAPI:
    """Simple KEGG REST client with basic argument validation."""

    BASE_URL = "https://rest.kegg.jp"
    _MIN_INTERVAL_SECONDS = 0.34

    def __init__(self) -> None:
        self.session = requests.Session()
        self._last_call_ts = 0.0

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_call_ts
        if elapsed < self._MIN_INTERVAL_SECONDS:
            time.sleep(self._MIN_INTERVAL_SECONDS - elapsed)
        self._last_call_ts = time.monotonic()

    def query(self, op, arg1, arg2=None, arg3=None):
        """Build and execute KEGG request."""
        args = [str(op), str(arg1)]
        if arg2 is not None:
            args.append(str(arg2))
        if arg3 is not None:
            args.append(str(arg3))

        request_url = f"{self.BASE_URL}/{'/'.join(args)}"
        self._throttle()
        results = self.session.get(request_url, timeout=60)
        results.raise_for_status()
        print(request_url)
        print(results.text)
        return results.text

    def kegg_info(self, database: str) -> str:
        return self.query("info", database)

    def kegg_list(self, database, org=None) -> str:
        if database in ("pathway", "module") and org:
            return self.query("list", database, org)
        if isinstance(database, list):
            if len(database) > 100:
                raise ValueError("Maximum number of databases is 100 for kegg list query")
            joined_database = "+".join(database)
            return self.query("list", joined_database)
        if org is not None:
            raise ValueError("Invalid database arg for kegg list request.")
        return self.query("list", database)

    def kegg_find(self, database, query, option=None) -> str:
        if isinstance(query, list):
            query = "+".join(query)
        if option is not None:
            if database in ("compound", "drug") and option in ("formula", "exact_mass", "mol_weight"):
                return self.query("find", database, query, option)
            raise ValueError("Invalid option arg for kegg find request.")
        return self.query("find", database, query)

    def kegg_get(self, dbentries, option=None) -> str:
        if isinstance(dbentries, list):
            if len(dbentries) > 10:
                raise ValueError("Maximum number of dbentries is 10 for kegg get query")
            dbentries = "+".join(dbentries)
        if option is not None:
            if option in ("aaseq", "ntseq", "mol", "kcf", "image", "kgml", "json"):
                return self.query("get", dbentries, option)
            raise ValueError("Invalid option arg for kegg get request.")
        return self.query("get", dbentries)

    def kegg_conv(self, target_db, source_db, option=None) -> str:
        if isinstance(source_db, list):
            source_db = "+".join(source_db)

        if option is not None and option not in ("turtle", "n-triple"):
            raise ValueError("Invalid option arg for kegg conv request.")

        gene_dbs = {"ncbi-gi", "ncbi-geneid", "uniprot"}
        chem_kegg_dbs = {"drug", "compound", "glycan"}
        chem_external_dbs = {"pubchem", "glycan"}
        condition1 = target_db in gene_dbs or source_db in gene_dbs
        condition2 = target_db in chem_kegg_dbs and source_db in chem_external_dbs
        condition3 = target_db in chem_external_dbs and source_db in chem_kegg_dbs
        if not (condition1 or condition2 or condition3):
            raise ValueError("Bad argument target_db or source_db for kegg conv request.")

        if option is not None:
            return self.query("conv", target_db, source_db, option)
        return self.query("conv", target_db, source_db)

    def kegg_link(self, target_db, source_db, option=None) -> str:
        if isinstance(source_db, list):
            source_db = "+".join(source_db)
        if option is not None and option not in ("turtle", "n-triple"):
            raise ValueError("Invalid option arg for kegg conv request.")
        if option is not None:
            return self.query("link", target_db, source_db, option)
        return self.query("link", target_db, source_db)
