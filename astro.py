import pyswisseph as se
from datetime import date
from historical_pattern_finder import PLANETS_HEBREW_MAP, SIGNS_HEBREW_MAP, find_historical_pattern

def calculate_houses(jd, lat, lon):
    """
    מחשב את בתי האסטרולוגיה עבור תאריך, קו רוחב וקו אורך נתונים.
    """
    house_system = b'P' # פלסידוס
    houses, ascmc = se.houses(jd, lat, lon, house_system)
    
    house_positions = {}
    for i, h_pos in enumerate(houses):
        house_positions[f"בית {i+1}"] = {
            "longitude": h_pos,
            "sign": SIGNS_HEBREW_MAP.get(int(h_pos // 30)),
            "degree_in_sign": h_pos % 30
        }
    return house_positions


def calculate_birth_chart(target_date, lat, lon):
    """
    מחשב מפת לידה מלאה: מיקומי כוכבים ומיקומי בתים.
    """
    se.set_ephe_path('swisseph_data/')
    se.swe_set_ephe_path('swisseph_data/')
    
    jd = se.julday(target_date.year, target_date.month, target_date.day)

    birth_chart = {}

    # חישוב מיקומי כוכבים
    for planet_he, planet_en in PLANETS_HEBREW_MAP.items():
        if planet_en != "thali":
            # עבור כוכבי לכת רגילים
            if planet_en == "north_node":
                flag = se.SEFLG_TRUEPOS
                lon, _, _ = se.swe_calc_ut(jd, se.NORT_NODE, flag)
            else:
                lon, _, _ = se.swe_calc_ut(jd, getattr(se, f"SE_{planet_en.upper()}"))
            
            sign_index = int(lon[0] // 30)
            sign_hebrew = SIGNS_HEBREW_MAP.get(sign_index)
            degree_in_sign = lon[0] % 30
            
            birth_chart[planet_he] = {
                "longitude": lon[0],
                "sign": sign_hebrew,
                "degree_in_sign": degree_in_sign
            }
    
    # הוספת תלי
    thali_data = find_historical_pattern("תלי", None, target_date, target_date)
    if thali_data:
        thali_pos = thali_data[0]
        birth_chart["תלי"] = {
            "longitude": thali_pos["longitude"],
            "sign": thali_pos["sign"],
            "degree_in_sign": thali_pos["degree_in_sign"]
        }

    # חישוב מיקומי בתים
    house_positions = calculate_houses(jd, lat, lon)
    birth_chart["בתים"] = house_positions
    
    return birth_chart
