import swisseph as se
from datetime import date, timedelta
import json

# מילון שמות כוכבים
PLANETS_HEBREW_MAP = {
    "שמש": se.SUN, "ירח": se.MOON, "מרקורי": se.MERCURY, "ונוס": se.VENUS,
    "מאדים": se.MARS, "יופיטר": se.JUPITER, "שבתאי": se.SATURN,
    "אורנוס": se.URANUS, "נפטון": se.NEPTUNE, "פלוטו": se.PLUTO,
    # הוספת ראש וזנב תלי
    "ראש תלי": se.TRUE_NODE,
    "זנב תלי": se.MEAN_NODE # ניתן להשתמש ב-MEAN_NODE או ב-TRUE_NODE
}

# מילון שמות מזלות בעברית
SIGNS_HEBREW_MAP = {
    0: "טלה", 1: "שור", 2: "תאומים", 3: "סרטן", 4: "אריה", 5: "בתולה",
    6: "מאזניים", 7: "עקרב", 8: "קשת", 9: "גדי", 10: "דלי", 11: "דגים"
}

# מילון היבטים עם מעלות ואורבים
ASPECTS_HEBREW_MAP = {
    "צמידות": {"degree": 0, "orb": 2},
    "היפוך": {"degree": 180, "orb": 2},
    "טרין": {"degree": 120, "orb": 2},
    "סקסטיל": {"degree": 60, "orb": 2},
    "ריבוע": {"degree": 90, "orb": 2}
}

def get_sign_id(longitude):
    """מחזירה את מספר המזל (0-11) על פי קו האורך."""
    return int(longitude / 30)

def find_historical_pattern(planet_name, sign_name, start_date, end_date, target_degree=None, degree_tolerance=0.5):
    """
    סורקת תקופה היסטורית ומחזירה את כל התאריכים שבהם כוכב מסוים נמצא במזל מסוים
    ובמעלה נתונים.
    """
    planet_id = PLANETS_HEBREW_MAP.get(planet_name)
    try:
        sign_id = [key for key, value in SIGNS_HEBREW_MAP.items() if value == sign_name][0]
    except (IndexError, KeyError):
        print(f"שם מזל לא חוקי: {sign_name}. אנא בחר מתוך הרשימה:")
        print(list(SIGNS_HEBREW_MAP.values()))
        return []

    if planet_id is None:
        print(f"שם כוכב לא חוקי: {planet_name}. אנא בחר מתוך הרשימה:")
        print(list(PLANETS_HEBREW_MAP.keys()))
        return []

    found_dates = []
    current_date = start_date
    delta = timedelta(days=1)
    
    degree_info = f" במעלה {target_degree}" if target_degree is not None else ""
    print(f"מתחיל סריקה: מתי {planet_name} היה במזל {sign_name}{degree_info} מ- {start_date} עד {end_date}...")
    
    while current_date <= end_date:
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        planet_pos, _ = se.calc_ut(jd_utc, planet_id)
        longitude = planet_pos[0]
        
        degree_in_sign = longitude % 30
        
        is_in_sign = get_sign_id(longitude) == sign_id
        is_in_degree = True
        if target_degree is not None:
            is_in_degree = abs(degree_in_sign - target_degree) <= degree_tolerance
        
        if is_in_sign and is_in_degree:
            found_dates.append({
                "date": str(current_date),
                "longitude": round(longitude, 2),
                "degree_in_sign": round(degree_in_sign, 2),
                "sign": sign_name
            })
            
        current_date += delta
        
    print("סריקה הסתיימה.")
    return found_dates

def find_historical_aspect(planet1_name, planet2_name, aspect_name, start_date, end_date):
    """
    סורקת תקופה היסטורית ומחזירה את כל התאריכים שבהם שני כוכבים נמצאים בהיבט מסוים.
    """
    planet1_id = PLANETS_HEBREW_MAP.get(planet1_name)
    planet2_id = PLANETS_HEBREW_MAP.get(planet2_name)
    aspect_info = ASPECTS_HEBREW_MAP.get(aspect_name)
    
    if planet1_id is None or planet2_id is None:
        print("שם כוכב לא חוקי.")
        return []
    
    if aspect_info is None:
        print("שם היבט לא חוקי. אנא בחר מתוך הרשימה:")
        print(list(ASPECTS_HEBREW_MAP.keys()))
        return []
        
    found_aspects = []
    current_date = start_date
    delta = timedelta(days=1)

    print(f"מתחיל סריקה: מתי {planet1_name} היה ב{aspect_name} ל{planet2_name} מ- {start_date} עד {end_date}...")

    while current_date <= end_date:
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        pos1, _ = se.calc_ut(jd_utc, planet1_id)
        pos2, _ = se.calc_ut(jd_utc, planet2_id)
        
        lon1 = pos1[0]
        lon2 = pos2[0]
        
        diff = abs(lon1 - lon2)
        if diff > 180:
            diff = 360 - diff
            
        if abs(diff - aspect_info["degree"]) <= aspect_info["orb"]:
            found_aspects.append({
                "date": str(current_date),
                "planet1": planet1_name,
                "planet2": planet2_name,
                "aspect": aspect_name,
                "orb": round(abs(diff - aspect_info["degree"]), 2)
            })
        
        current_date += delta

    print("סריקה הסתיימה.")
    return found_aspects

# --- קוד לבדיקת הפונקציה ---
if __name__ == "__main__":
    # הגדרת טווח קצר לבדיקות מהירות
    start_date_short = date(2023, 1, 1)
    end_date_short = date(2025, 1, 1)

    # בדיקה 1: מתי ראש תלי היה במזל טלה?
    print("--- חיפוש: ראש תלי במזל טלה ---")
    node_in_aries = find_historical_pattern("ראש תלי", "טלה", start_date_short, end_date_short)
    if node_in_aries:
        print(f"נמצאו {len(node_in_aries)} אירועים.")
        print(json.dumps(node_in_aries[:5], indent=2, ensure_ascii=False) + "...")
    else:
        print("לא נמצאו אירועים בטווח התאריכים הקצר.")

    print("\n" + "="*50 + "\n")
    
    # בדיקה 2: מתי יופיטר היה בצמידות לשבתאי? (התרחש בדצמבר 2020)
    print("--- חיפוש: צמידות יופיטר-שבתאי בטווח קצר ---")
    jupiter_saturn_conjunctions = find_historical_aspect("יופיטר", "שבתאי", "צמידות", start_date_short, end_date_short)
    if jupiter_saturn_conjunctions:
        print(f"נמצאו {len(jupiter_saturn_conjunctions)} אירועים.")
        print(json.dumps(jupiter_saturn_conjunctions, indent=2, ensure_ascii=False))
    else:
        print("לא נמצאו אירועים בטווח התאריכים הקצר.")