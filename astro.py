import pyswisseph as se
from datetime import date, datetime
from typing import Dict, Any, List, Optional
import itertools

# הגדרת נתיב לקבצי האפמריס (חובה עבור swisseph)
se.set_ephe_path('ephe')

# מילון שמות כוכבים
PLANETS_HEBREW_MAP = {
    "שמש": se.SUN, "ירח": se.MOON, "מרקורי": se.MERCURY, "ונוס": se.VENUS,
    "מאדים": se.MARS, "יופיטר": se.JUPITER, "שבתאי": se.SATURN,
    "אורנוס": se.URANUS, "נפטון": se.NEPTUNE, "פלוטו": se.PLUTO,
    "ראש תלי": se.TRUE_NODE, "זנב תלי": se.MEAN_NODE
}

# מילון מזלות עם טווח המעלות שלהם
SIGNS_HEBREW_MAP_DEGREE_RANGES = {
    "טלה": (0, 30), "שור": (30, 60), "תאומים": (60, 90),
    "סרטן": (90, 120), "אריה": (120, 150), "בתולה": (150, 180),
    "מאזניים": (180, 210), "עקרב": (210, 240), "קשת": (240, 270),
    "גדי": (270, 300), "דלי": (300, 330), "דגים": (330, 360)
}

# רשימת שמות המזלות בסדר (לשימוש בפונקציות אחרות שאולי תלויות בסדר)
SIGNS_HEBREW_LIST = [
    "טלה", "שור", "תאומים", "סרטן", "אריה", "בתולה",
    "מאזניים", "עקרב", "קשת", "גדי", "דלי", "דגים"
]

# הגדרת אורב כללי לתבניות גאומטריות (חשוב שיהיה כאן, כי זהו מודול הליבה לחישובים)
GEOMETRIC_PATTERN_ORB = 5 # מעלות

# --- פונקציות עזר חדשות/מעודכנות לטיפול בתאריכים ---

def get_date_from_components(year: int, month: int, day: int) -> datetime:
    """
    יוצר אובייקט datetime מתאריך הנתון.
    פונקציה זו מקבלת שנה בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    לדוגמה: שנת 1 לפנה"ס היא 0, שנת 2 לפנה"ס היא -1.
    """
    # datetime.date ו-datetime.datetime מטפלים בשנים שליליות כפי שמצופה עבור שנים אסטרונומיות.
    # עם זאת, עבור 0000 יש בעיה, לכן נמיר אם שנה היא 0.
    # Swiss Ephemeris (se.utc_to_jd) מצפה לשנה 0 עבור 1 לפנה"ס.
    # Python datetime לא תומך בשנת 0 (יראה אותה כ-1 לספירה כנראה, או יזרוק שגיאה תלוי בגרסה).
    # לכן, נדאג להעביר ל-se.utc_to_jd את השנה כפי שקיבלנו (שלילי/חיובי)
    # וליצור אובייקט datetime רק במידת הצורך (לדוגמה, עבור שמירה או הצגה).
    # עבור ה-datetime, נמנע מ-0000 אם יגיע
    if year == 0: # 1 BCE is year 0 in astronomical numbering, but datetime doesn't like 0
        # For internal datetime object, treat year 0 as 1 AD and adjust the flag if needed,
        # but for swisseph, pass the astronomical year directly.
        # This function is primarily for internal datetime object creation logic consistency.
        # The frontend should already send astronomical year for BCE.
        # If year 0 comes from frontend (representing 1 BCE), it should be handled internally.
        # For direct datetime creation, a dummy non-zero year is sometimes used with careful handling.
        # However, for swisseph, 0 is valid.
        pass # The logic below will use the year directly.

    try:
        return datetime(year, month, day)
    except ValueError as e:
        # Handle cases like invalid day for month etc.
        raise ValueError(f"Invalid date components: Year={year}, Month={month}, Day={day}. Error: {e}")


# --- פונקציות אסטרולוגיות ליבה ---

def get_sign_from_longitude(longitude: float) -> str:
    """מחזירה את שם המזל לפי קו אורך (מ-0 עד 360 מעלות)."""
    for sign, (start_deg, end_deg) in SIGNS_HEBREW_MAP_DEGREE_RANGES.items():
        if start_deg <= longitude < end_deg:
            return sign
    return "לא ידוע"

def get_angular_distance(lon1: float, lon2: float) -> float:
    """
    מחשבת את ההפרש הזוויתי הקצר ביותר בין שני קווי אורך.
    """
    diff = abs(lon1 - lon2)
    if diff > 180:
        diff = 360 - diff
    return diff

def find_all_grand_trines(positions: Dict[str, float]) -> List[Dict[str, Any]]:
    """
    מזהה את כל משולשי היסוד הגדולים (Grand Trines) מתוך מיקומי כוכבים.
    Grand Trine: 3 כוכבים, כל אחד במרחק של כ-120 מעלות מהאחרים.
    """
    planets_with_positions = [(p_name, p_lon) for p_name, p_lon in positions.items() if p_lon is not None]
    found_trines = []
    
    # בודק קומבינציות של 3 כוכבים
    for p1, p2, p3 in itertools.combinations(planets_with_positions, 3):
        p1_name, p1_lon = p1
        p2_name, p2_lon = p2
        p3_name, p3_lon = p3

        d12 = get_angular_distance(p1_lon, p2_lon)
        d23 = get_angular_distance(p2_lon, p3_lon)
        d31 = get_angular_distance(p3_lon, p1_lon)

        is_trine_12 = abs(d12 - 120) <= GEOMETRIC_PATTERN_ORB
        is_trine_23 = abs(d23 - 120) <= GEOMETRIC_PATTERN_ORB
        is_trine_31 = abs(d31 - 120) <= GEOMETRIC_PATTERN_ORB

        if is_trine_12 and is_trine_23 and is_trine_31:
            found_trines.append({
                "name": "Grand Trine",
                "planets": sorted([p1_name, p2_name, p3_name]),
                "degrees": {
                    f"{p1_name}-{p2_name}": round(d12, 2),
                    f"{p2_name}-{p3_name}": round(d23, 2),
                    f"{p3_name}-{p1_name}": round(d31, 2)
                }
            })
    return found_trines


def find_star_of_david(positions: Dict[str, float]) -> Optional[Dict[str, Any]]:
    """
    מזהה תבנית "מגן דוד" (Star of David) מתוך מיקומי כוכבים.
    מגן דוד: שני משולשי יסוד גדולים (Grand Trines) שיוצרים ביניהם 3 אופוזיציות.
    דורש לפחות 6 כוכבים.
    """
    planets_with_positions = [(p_name, p_lon) for p_name, p_lon in positions.items() if p_lon is not None]

    if len(planets_with_positions) < 6:
        return None

    all_trines = find_all_grand_trines(positions)
    
    # חפש זוגות של Grand Trines
    for gt1, gt2 in itertools.combinations(all_trines, 2):
        planets_gt1 = set(gt1["planets"])
        planets_gt2 = set(gt2["planets"])
        
        # וודא שהכוכבים בשני הטרינים שונים לחלוטין ושסך הכוכבים הוא 6
        if len(planets_gt1.union(planets_gt2)) != 6 or len(planets_gt1.intersection(planets_gt2)) > 0:
            continue

        sod_planets = list(planets_gt1.union(planets_gt2))
        
        current_sod_positions = {p_name: positions[p_name] for p_name in sod_planets}
        current_sod_planets_list = [(p, current_sod_positions[p]) for p in sod_planets]

        found_oppositions = []
        # נבדוק אופוזיציות בין כוכבים מהטרין הראשון לכוכבים מהטרין השני
        for p_a_name in planets_gt1:
            for p_b_name in planets_gt2:
                dist = get_angular_distance(positions[p_a_name], positions[p_b_name])
                if abs(dist - 180) <= GEOMETRIC_PATTERN_ORB:
                    found_oppositions.append(tuple(sorted([p_a_name, p_b_name])))

        # מגן דוד דורש בדיוק 3 אופוזיציות ייחודיות (כל כוכב משתתף באופוזיציה אחת)
        # נסנן כפילויות של אופוזיציות (כיוון A-B זהה ל-B-A)
        unique_found_oppositions = list(set(found_oppositions))

        # מגן דוד דורש שכל 6 הכוכבים יהיו מעורבים באופוזיציות, ושיהיו 3 אופוזיציות בדיוק.
        # בנוסף, כל כוכב חייב להופיע באופוזיציה אחת בלבד.
        all_opposition_planets = set()
        for op in unique_found_oppositions:
            all_opposition_planets.add(op[0])
            all_opposition_planets.add(op[1])

        if len(unique_found_oppositions) == 3 and len(all_opposition_planets) == 6:
            return {
                "name": "Star of David",
                "planets": sorted(sod_planets),
                "trines": [gt1, gt2],
                "oppositions": unique_found_oppositions
            }
            
    return None

def _calculate_houses(jd: float, lat: float, long_geo: float) -> Dict[str, Any]:
    """
    מחשב את בתי האסטרולוגיה עבור תאריך, קו רוחב וקו אורך נתונים.
    """
    house_system = b'P' # פלסידוס (Placidus)
    houses, ascmc = se.houses(jd, lat, long_geo, house_system)

    house_positions = {}
    for i, h_pos in enumerate(houses):
        longitude = h_pos % 360
        sign = get_sign_from_longitude(longitude)
        degree_in_sign = int(longitude) % 30

        house_positions[f"בית {i+1}"] = {
            "longitude": round(longitude, 2),
            "sign": sign,
            "degree_in_sign": degree_in_sign
        }
    
    house_positions["אופק (Ascendant)"] = {
        "longitude": round(ascmc[0] % 360, 2),
        "sign": get_sign_from_longitude(ascmc[0] % 360),
        "degree_in_sign": int(ascmc[0] % 360) % 30
    }
    house_positions["רום שמיים (Midheaven)"] = {
        "longitude": round(ascmc[1] % 360, 2),
        "sign": get_sign_from_longitude(ascmc[1] % 360),
        "degree_in_sign": int(ascmc[1] % 360) % 30
    }

    return house_positions

def calculate_birth_chart(year: int, month: int, day: int, lat: float, long_geo: float) -> Dict[str, Any]:
    """
    מחשב מפת לידה מלאה: מיקומי כוכבים, מיקומי בתים, ותבניות גאומטריות.
    מקבל שנה בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    מחזירה מילון מובנה הכולל את כל הנתונים.
    **מעודכן**: כולל מהירויות כוכבים בפלט.
    """
    # se.utc_to_jd מקבל שנה כפי שהיא, כולל שליליות עבור שנים אסטרונומיות.
    # שנה 0 היא שנת 1 לפנה"ס. שנה -1 היא שנת 2 לפנה"ס וכן הלאה.
    jd_utc = se.utc_to_jd(year, month, day, 0, 0, 0)[1]

    # יצירת אובייקט תאריך עבור שמירה/הצגה (ייתכן ותהיה פחות עמידה לשנים מאוד שליליות ב-Python מובנה)
    # לצורך הצגה, נשתמש בפורמט שיתאים יותר להצגה ב-frontend (לדוגמה, 1563 לפנה"ס)
    # ולא בשנה אסטרונומית (לדוגמה, -1562).
    # הפרונטאנד יטפל בהצגה הנכונה של BCE, אז כאן נשמור את השנה כמו שהיא,
    # או שנה סטנדרטית לצורך הדפסה בלבד בבדיקות עצמיות.
    display_date_string = f"{day}/{month}/{year}" # השנה כאן היא כבר אסטרונומית
    if year <= 0: # אם זו שנה אסטרונומית 0 או שלילית
        display_year_bce = -(year - 1) # המרה לשנת לפנה"ס רגילה
        display_date_string = f"{day}/{month}/{display_year_bce} לפנה'ס"

    chart_data = {
        "date": display_date_string, # נשמור את התאריך בפורמט קריא יותר (גם אם השנה אסטרונומית)
        "latitude": lat,
        "longitude": long_geo,
        "planet_positions": {},
        "house_positions": {},
        "detected_patterns": [] # רשימה חדשה לתבניות מזוהות
    }

    raw_planet_longitudes = {} # נשמור קווי אורך כדי להעביר לפונקציות זיהוי תבניות
    for planet_name_hebrew, planet_id in PLANETS_HEBREW_MAP.items():
        try:
            # שינוי: קוראים עם se.SWE_FLG_SPEED כדי לקבל גם מהירות
            position, _ = se.calc_ut(jd_utc, planet_id, se.SWE_FLG_SPEED)
            
            longitude = position[0] % 360
            speed = position[3] # מהירות הכוכב במעלות ליום
            
            raw_planet_longitudes[planet_name_hebrew] = longitude # שמירת קו האורך הגולמי

            sign = get_sign_from_longitude(longitude)
            degree_in_sign = int(longitude) % 30

            chart_data["planet_positions"][planet_name_hebrew] = {
                "longitude": round(longitude, 2),
                "sign": sign,
                "degree_in_sign": degree_in_sign,
                "speed": round(speed, 4), # הוספה של מהירות הכוכב
                "is_retrograde": speed < 0 # סימון נסיגה אם המהירות שלילית
            }
        except Exception as e:
            # print(f"שגיאה בחישוב מיקום {planet_name_hebrew} בתאריך {year}/{month}/{day}: {e}") # לוגינג של שגיאות
            chart_data["planet_positions"][planet_name_hebrew] = {
                "error": f"שגיאה בחישוב מיקום: {e}"
            }
    
    # חישוב מיקומי בתים
    chart_data["house_positions"] = _calculate_houses(jd_utc, lat, long_geo)

    # זיהוי תבניות גאומטריות מתוך מיקומי הכוכבים שחושבו
    if len(raw_planet_longitudes) >= 3: # Grand Trine דורש לפחות 3 כוכבים
        found_gts = find_all_grand_trines(raw_planet_longitudes)
        if found_gts:
            chart_data["detected_patterns"].extend(found_gts)
    
    if len(raw_planet_longitudes) >= 6: # מגן דוד דורש לפחות 6 כוכבים
        sod = find_star_of_david(raw_planet_longitudes)
        if sod:
            chart_data["detected_patterns"].append(sod)

    return chart_data

# פונקציות dummy לשלמות, ייתכן ותצטרך לממש אותן במידה וצריך
def get_planets_and_signs_data():
    return {} # למלא בנתונים אם נדרש

def get_aspects_data():
    return {} # למלא בנתונים אם נדרש

def calculate_moon_phase(year: int, month: int, day: int) -> Dict[str, Any]:
    """
    מחשב את שלב הירח לתאריך נתון.
    (זוהי פונקציית דמה, יש לממש את הלוגיקה בפועל)
    """
    # כאן תבוא לוגיקת חישוב שלב הירח באמצעות swisseph
    # לדוגמה, מציאת ההפרש בין השמש לירח
    return {"phase": "ירח מלא", "illumination": 99.5}


# בדיקה עצמית של הפונקציה (רק אם מריצים את הקובץ ישירות)
if __name__ == "__main__":
    import json
    
    test_lat = 31.7683 # ירושלים
    test_long = 35.2137 # ירושלים

    print("--- חישוב מפת לידה לדוגמה (ירושלים, 1986-04-13) ---")
    chart = calculate_birth_chart(1986, 4, 13, test_lat, test_long)
    print(json.dumps(chart, indent=2, ensure_ascii=False))

    print("\n--- חישוב מפת לידה לדוגמה (17/08/1987) עם בדיקת תבניות ---")
    # הגדל אורב לצורך הדגמה (בד"כ אורבים קטנים יותר)
    original_orb = GEOMETRIC_PATTERN_ORB
    GEOMETRIC_PATTERN_ORB = 8 
    chart_with_patterns = calculate_birth_chart(1987, 8, 17, test_lat, test_long)
    print(json.dumps(chart_with_patterns, indent=2, ensure_ascii=False))
    GEOMETRIC_PATTERN_ORB = original_orb # החזרת האורב

    print("\n--- חישוב מפת לידה לדוגמה לשנה לפני הספירה (1563 לפנה''ס, 01/07) ---")
    # 1563 לפנה"ס היא שנה -1562 בפורמט אסטרונומי
    chart_bce = calculate_birth_chart(-1562, 7, 1, test_lat, test_long)
    print(json.dumps(chart_bce, indent=2, ensure_ascii=False))

    print("\n--- חישוב מפת לידה לדוגמה לשנה 1 לפנה''ס (שנה 0 אסטרונומית) ---")
    # 1 לפנה"ס היא שנה 0 בפורמט אסטרונומי
    chart_1bce = calculate_birth_chart(0, 1, 1, test_lat, test_long)
    print(json.dumps(chart_1bce, indent=2, ensure_ascii=False))
