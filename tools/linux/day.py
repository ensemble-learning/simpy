from pyowm import OWM

API_key = 'bc892b51e40a9bd87478010fd7430b4d'
owm = OWM(API_key)
obs = owm.weather_at_coords(34.15, -118.14)
w = obs.get_weather()
t = w.get_temperature(unit='celsius')
t_max = t['temp_max']
t_min = t['temp_min']
print t_max, t_min


