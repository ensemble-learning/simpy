import time

Second2Day = 60*60*24
items = []

today = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 
today = today[:10]
t0 = time.mktime(time.strptime(today, "%Y-%m-%d"))

day = "2016-10-10"
title = "JACS revision"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

day = "2016-10-20"
title = "PCCP revision"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

day = "2016-10-17"
title = "wedding anniversary"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

day = "2016-11-19"
title = "PNAS revision"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

day = "2016-09-20"
title = "JACS submission (Mark)"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

day = "2016-08-25"
title = "Science submission (Li)"

t1 = time.mktime(time.strptime(day, "%Y-%m-%d"))
n = int((t1 - t0)/Second2Day)
items.append([n, title])

return items


