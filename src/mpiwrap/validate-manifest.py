#!/usr/bin/env python3
"""Check that the MPI support matrix and marshalling metadata stay in lockstep."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MATRIX = ROOT / "support-matrix.tsv"
METADATA = ROOT / "argument-metadata.json"
EXCLUDED = {"excluded-mpi5-removed", "excluded-deprecated", "excluded-non-C-binding"}


def fail(message):
    raise SystemExit("manifest error: " + message)


def main():
    with MATRIX.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    functions = {
        row["symbol"]: row for row in rows
        if row["kind"] == "function" and row["symbol"].startswith("MPI_")
        and row["in_scope"] not in EXCLUDED
        and row["symbol"] != "MPI_Comm_from_native"
    }
    metadata = json.loads(METADATA.read_text())
    if metadata.get("schema_version") != 3:
        fail("unsupported argument-metadata schema")
    entries = metadata.get("functions", [])
    names = [entry.get("name") for entry in entries]
    if len(names) != len(set(names)):
        fail("duplicate function metadata names")
    missing = sorted(set(functions) - set(names))
    extra = sorted(set(names) - set(functions))
    if missing or extra:
        fail(f"function coverage mismatch; missing={missing}, extra={extra}")
    count = 0
    for entry in entries:
        name = entry["name"]
        row = functions[name]
        if entry.get("signature") != row["api_signature_or_representation"]:
            fail(f"signature differs between matrix and argument metadata: {name}")
        if not entry.get("provider_behavior") or not entry.get("test_id"):
            fail(f"provider behavior and semantic test are required: {name}")
        routes = entry.get("provider_routes", {})
        for provider in ("native_real_mpi", "internal_thread_mpi", "injected_real_mpi"):
            if not routes.get(provider, {}).get("status"):
                fail(f"provider route status is required: {name}/{provider}")
        if routes["internal_thread_mpi"]["status"] == "unavailable" and not routes["internal_thread_mpi"].get("reason"):
            fail(f"unavailable provider operation needs a reason: {name}")
        if routes["internal_thread_mpi"]["status"] == "narrow-collective-emulation":
            if not routes["internal_thread_mpi"].get("prerequisites") or not routes["internal_thread_mpi"].get("failure"):
                fail(f"emulation must list prerequisites and failure behavior: {name}")
        args = entry.get("arguments", [])
        if len(args) != len({arg.get("name") for arg in args}):
            fail(f"duplicate argument metadata: {name}")
        for arg in args:
            pointer = "*" in arg["c_type"] or "[" in arg["c_type"] or arg["c_type"] == "..."
            if arg.get("direction") not in {"in", "out", "inout", "variadic"}:
                fail(f"invalid direction {name}.{arg.get('name')}: {arg.get('direction')}")
            if not arg.get("conversion"):
                fail(f"missing conversion: {name}.{arg['name']}")
            if pointer and not arg.get("extent"):
                fail(f"pointer argument requires an explicit extent/lifetime: {name}.{arg['name']}")
            if arg.get("conversion") in {"array", "handle_array", "status_array", "integer_array"}:
                if arg["direction"] in {"out", "inout"} and not arg.get("valid_on_return"):
                    fail(f"output array must specify valid entries on return: {name}.{arg['name']}")
            if arg.get("conversion") == "typed_callback":
                matching = [callback for callback in entry.get("callbacks", [])
                            if callback.get("argument") == arg["name"]]
                if len(matching) != 1:
                    fail(f"typed callback requires exactly one lifetime/bridge record: {name}.{arg['name']}")
                for key in ("native_signature", "registration", "required_bridge", "provider_mode"):
                    if not matching[0].get(key):
                        fail(f"callback metadata missing {key}: {name}.{arg['name']}")
            count += 1
    constants = json.loads((ROOT / "constant-metadata.json").read_text())
    if constants.get("schema_version") != 1:
        fail("unsupported constant-metadata schema")
    required_constants = {row["symbol"] for row in rows
                          if row["kind"] == "constant" and row["in_scope"] not in EXCLUDED}
    constant_entries = constants.get("constants", [])
    constant_names = [item.get("standard_api_name") for item in constant_entries]
    if set(constant_names) != required_constants or len(constant_names) != len(set(constant_names)):
        fail("constant coverage mismatch or duplicates")
    for item in constant_entries:
        if item.get("class") not in {"provider_handle", "provider_integer", "compile_time_integer",
                                      "buffer_sentinel", "status_sentinel", "status_array_sentinel"}:
            fail(f"constant needs explicit value class: {item.get('standard_api_name')}")
        if not item.get("translation") or not item.get("thread_mpi_2026_4", {}).get("status"):
            fail(f"constant needs a conversion and tMPI support status: {item.get('standard_api_name')}")
    type_data = json.loads((ROOT / "type-metadata.json").read_text())
    if type_data.get("schema_version") != 1:
        fail("unsupported type-metadata schema")
    required_types = {row["symbol"] for row in rows if row["kind"] == "type"}
    type_entries = type_data.get("types", [])
    type_names = [item.get("name") for item in type_entries]
    if set(type_names) != required_types or len(type_names) != len(set(type_names)):
        fail("type coverage mismatch or duplicates")
    for item in type_entries:
        if item.get("category") == "unclassified":
            fail(f"type needs explicit conversion/storage classification: {item.get('name')}")
        if item.get("category") == "callback-type" and not all(
            item.get(key) for key in ("signature", "bridge", "provider")
        ):
            fail(f"callback type needs exact signature, bridge, and provider policy: {item.get('name')}")
    print(f"validated {len(entries)} function rows/{count} arguments, {len(constant_entries)} constants, and {len(type_entries)} types")


if __name__ == "__main__":
    main()
