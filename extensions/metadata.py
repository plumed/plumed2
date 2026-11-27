"""Source declarations shared by module inventory and documentation tools."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/maketools"))
from module_header import check_obsolete, header_field, read_header

# check_obsolete and header_field are re-exported for the inventory and the manual
__all__ = ["check_obsolete", "header_field", "read_header", "read_metadata"]


def read_metadata(path):
    """Load YAML without duplicate keys, custom tags, aliases or merge keys."""
    import yaml

    class Loader(yaml.SafeLoader):
        def compose_node(self, parent, index):
            if self.check_event(yaml.AliasEvent):
                raise ValueError(f"{path}: YAML aliases are not supported")
            return super().compose_node(parent, index)

        def construct_mapping(self, node, deep=False):
            result = {}
            for key_node, value_node in node.value:
                key = self.construct_object(key_node, deep=deep)
                if not isinstance(key, str) or key in result:
                    raise ValueError(f"{path}: non-string or duplicate key {key!r}")
                result[key] = self.construct_object(value_node, deep=deep)
            return result

    text = Path(path).read_text()
    try:
        data = yaml.load(text, Loader=Loader)
    except yaml.YAMLError as error:
        raise ValueError(f"{path}: {error}") from error
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    if Path(path).name == "module.yml":
        header = read_header(path, text)
        if any(data.get(key) != value for key, value in header.items()):
            raise ValueError(f"{path}: YAML and native header disagree")
    return data
