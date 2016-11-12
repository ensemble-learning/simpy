import os
from pyowm import OWM
import datetime
import time

class Weather():
    def __init__(self,):
        self.city = ""
        self.status = ""
        self.temp_min = 0.0
        self.temp_max = 0.0

class Todos():
    def __init__(self,):
        self.items = []

class Events():
    def __init__(self,):
        self.items = []
        self.items_clean = []
    def parse(self,):
        n = 0
        for i in self.items:
            if n > 0:
                tokens = i.strip()
                date = tokens[:10]
                time = tokens[11:16]
                title = tokens[26:]
                item = date + " " + time + " " + title + "\n"
                self.items_clean.append(item)
            n += 1

def get_weather(weather, today):
    # get the weather
    API_key = 'bc892b51e40a9bd87478010fd7430b4d'
    owm = OWM(API_key)
    obs = owm.weather_at_coords(34.15, -118.14)
    w = obs.get_weather()
    weather.city = obs.get_location()
    weather.status = w.get_status() # status
    t = w.get_temperature(unit='celsius')
    weather.temp_max = t['temp_max']
    weather.temp_min= t['temp_min']

def get_todos(t):
    f = os.popen("/home/tao/anaconda2/bin/habitica todos")
    for i in f:
        t.items.append(i)
    f.close()

def get_events(e):
    os.chdir("/home/tao/Soft/google/calendar")
    f = os.popen("/home/tao/anaconda2/bin/python quickstart.py")
    for i in f:
        e.items.append(i)
    f.close()
    e.parse()

def get_days_matter(items):

    Second2Day = 60*60*24

    today = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 
    today = today[:10]
    t0 = time.mktime(time.strptime(today, "%Y-%m-%d"))

    day = "2016-11-04"
    title = "J1 no objective letter"

    t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
    n = int((t1 - t0)/Second2Day)
    if n > 0:
        line = "%s is due in %d days\n"%(title, n)
    if n == 0:
        line = "%s is today\n"%(title)
    if n < 0:
        line = "%s has been %d days\n"%(title, -n)
    items.append(line)

    day = "2016-10-20"
    title = "PCCP revision"

    t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
    n = int((t1 - t0)/Second2Day)
    if n > 0:
        line = "%s is due in %d days\n"%(title, n)
    if n == 0:
        line = "%s is today\n"%(title)
    if n < 0:
        line = "%s has been %d days\n"%(title, -n)
    items.append(line)

    day = "2017-10-17"
    title = "wedding anniversary"

    t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
    n = int((t1 - t0)/Second2Day)
    if n > 0:
        line = "%s is due in %d days\n"%(title, n)
    if n == 0:
        line = "%s is today\n"%(title)
    if n < 0:
        line = "%s has been %d days\n"%(title, -n)
    items.append(line)

    day = "2016-11-19"
    title = "PNAS revision"

    t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
    n = int((t1 - t0)/Second2Day)
    if n > 0:
        line = "%s is due in %d days\n"%(title, n)
    if n == 0:
        line = "%s is today\n"%(title)
    if n < 0:
        line = "%s has been %d days\n"%(title, -n)
    items.append(line)

    day = "2016-09-20"
    title = "JACS submission (Mark)"

    t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
    n = int((t1 - t0)/Second2Day)
    if n > 0:
        line = "%s is due in %d days\n"%(title, n)
    if n == 0:
        line = "%s is today\n"%(title)
    if n < 0:
        line = "%s has been %d days\n"%(title, -n)
    items.append(line)

    return items

def output_daily(today, weather, todos, events, days_matter):
    o = open("tcheng-%s.txt"%today, "w")
    o.write("%s\n"%('-'*32))
    o.write("%s\n"%today)
    o.write("Weather: %s "%weather.status)
    o.write("from %.2f to %.2f\n"%(weather.temp_min, weather.temp_max))
    o.write("%s\n"%('-'*32))
    o.write("My todos:\n")
    o.write("%s\n"%('-'*32))
    for i in todos.items:
        o.write(i)
    o.write("My current events:\n")
    o.write("%s\n"%('-'*32))
    for i in range(6):
        o.write(events.items_clean[i])
    o.write("Days matter:\n")
    o.write("%s\n"%('-'*32))
    for i in days_matter:
        o.write(i)
    o.close()

    time.sleep(1)
    os.system('mail -s "tcheng-%s" chengtaoliterature@gmail.com < tcheng-%s.txt'%
            (today, today))

def main():
    today = datetime.date.today()
    weather = Weather()
    todos = Todos()
    events = Events()
    days_matter = []

    #get_weather(weather, today)
    get_todos(todos)
    get_events(events)
    get_days_matter(days_matter)
    output_daily(today, weather, todos, events, days_matter)
    
main()

