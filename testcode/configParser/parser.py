import ConfigParser

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

Config = ConfigParser.ConfigParser()
Config.read("config")
sections = Config.sections()

for i in sections:
    dict = ConfigSectionMap(i)
    print dict

cfgfile = open("next.ini",'w')

# add the settings to the structure of the file, and lets write it out...
Config.add_section('Person')
Config.set('Person','HasEyes',True)
Config.set('Person','Age', 50)
Config.write(cfgfile)
cfgfile.close()
