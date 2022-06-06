import json
import sys
import os
import re

if len(sys.argv) != 2:
    print("python %s configue_file" % sys.argv[0])
    exit()


def resolveConfig(configure_file):
    config_dict = {}
    with open(configure_file, 'r') as inopen:
        for line in inopen:
            line = line.strip()
            if line.startswith("//"):
                items = line[2:]
                config_dict.update({items: {}})
            elif line != "" and not line.startswith('#'):
                parameter = line.split(" = ")[0]
                value = line.split(" = ")[1].replace('None', '')
                try:
                    basic = re.findall(r'\$\{(.+)\}', value)[0]
                    value = value.replace(
                        '${%s}' % basic, config_dict.get('basic').get(basic))
                except IndexError:
                    pass
                config_dict.get(items).update({parameter: value})
    with open("%s.json" % configure_file, "w") as otopen:
        otopen.write(json.dumps(config_dict))


resolveConfig(sys.argv[1])
