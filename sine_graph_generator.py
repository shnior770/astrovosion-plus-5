import swisseph as se
from datetime import timedelta
from historical_pattern_finder import PLANETS_HEBREW_MAP, ASPECTS_HEBREW_MAP

def generate_sine_graph_data(planet_names, start_date, end_date, aspects_to_find=None):
    """
    מייצר נתוני מיקום כוכבים עבור גרף סינוסי, כולל איתור היבטים.
    """
    se.set_ephe_path('swisseph_data/')
    se.swe_set_ephe_path('swisseph_data/')
    
    graph_data = []
    
    # המרת שמות הכוכבים לאנגלית לצורך swisseph
    planets_en = {p_he: p_en for p_en, p_he in PLANETS_HEBREW_MAP.items()}
    
    current_date = start_date
    while current_date <= end_date:
        jd = se.julday(current_date.year, current_date.month, current_date.day)
        
        daily_data = {"date": str(current_date), "positions": {}, "aspects": []}
        
        # חישוב מיקום הכוכבים
        planet_positions = {}
        for planet_he in planet_names:
            planet_en = planets_en.get(planet_he)
            if planet_en:
                if planet_en == "north_node":
                    lon, _, _ = se.swe_calc_ut(jd, se.NORT_NODE, se.SEFLG_TRUEPOS)
                else:
                    lon, _, _ = se.swe_calc_ut(jd, getattr(se, f"SE_{planet_en.upper()}"))
                
                planet_positions[planet_he] = lon[0]
                daily_data["positions"][planet_he] = lon[0]
        
        # חישוב היבטים אם נדרש
        if aspects_to_find and len(planet_names) >= 2:
            found_aspects = []
            for i in range(len(planet_names)):
                for j in range(i + 1, len(planet_names)):
                    planet1_name = planet_names[i]
                    planet2_name = planet_names[j]
                    lon1 = planet_positions.get(planet1_name)
                    lon2 = planet_positions.get(planet2_name)
                    
                    if lon1 is not None and lon2 is not None:
                        for aspect_name, aspect_data in ASPECTS_HEBREW_MAP.items():
                            if aspect_name in aspects_to_find:
                                orb = 1  # אורב של מעלה אחת
                                diff = abs(lon1 - lon2)
                                if diff > 180:
                                    diff = 360 - diff
                                
                                if abs(diff - aspect_data["degrees"]) <= orb:
                                    found_aspects.append({
                                        "planets": f"{planet1_name}-{planet2_name}",
                                        "aspect": aspect_name,
                                        "degrees": aspect_data["degrees"],
                                        "orb": abs(diff - aspect_data["degrees"])
                                    })
            daily_data["aspects"] = found_aspects

        graph_data.append(daily_data)
        current_date += timedelta(days=1)
        
    return graph_data