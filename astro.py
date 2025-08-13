import pyswisseph as se
from datetime import date, datetime
from typing import Dict, Any, List, Optional
import itertools
import os # **חשוב**: ייבוא מודול os

# הגדרת נתיב לקבצי האפמריס (חובה עבור swisseph)
# **חשוב**: נשתמש במשתנה סביבה כדי לקבל נתיב מוחלט בתוך הקונטיינר
se.set_ephe_path(os.environ.get('SWISSEPH_PATH_ABS', 'ephe'))
# אם משתנה הסביבה 'SWISSEPH_PATH_ABS' לא מוגדר (לדוגמה, בהרצה מקומית ללא דוקר),
# הוא יחזור לברירת המחדל 'ephe'.

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

# --- חדש: הגדרת היבטים עם מעלות ואורבים (לפי דרישת המפתח המקורית) ---
ASPECTS_DEGREE_ORB_CONFIG = {
    "צמידות": {"degree": 0, "orb": 8, "type": "major"},
    "היפוך": {"degree": 180, "orb": 8, "type": "major"},
    "טרין": {"degree": 120, "orb": 8, "type": "major"},
    "ריבוע": {"degree": 90, "orb": 8, "type": "major"},
    "סקסטיל": {"degree": 60, "orb": 6, "type": "soft"},
    "קווינקונקס": {"degree": 150, "orb": 3, "type": "minor"} # נוסף היבט Quincunx
}

# --- פונקציות עזר חדשות/מעודכנות לטיפול בתאריכים ---

def get_date_from_components(year: int, month: int, day: int) -> datetime:
    """
    יוצר אובייקט datetime מתאריך הנתון.
    פונקציה זו מקבלת שנה בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    לדוגמה: שנת 1 לפנה"ס היא 0, שנת 2 לפנה"ס היא -1.
    """
    if year == 0:
        pass 
    try:
        return datetime(year, month, day)
    except ValueError as e:
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

# --- חדש: פונקציה לחישוב כל ההיבטים בין זוגות כוכבים ---
def calculate_all_aspects(
    planet_positions_raw: Dict[str, float]
) -> List[Dict[str, Any]]:
    """
    מחשבת את כל ההיבטים המשמעותיים בין כל זוגות כוכבי הלכת הנתונים.
    משתמשת בהגדרות ASPECT_DEGREES_ORB_CONFIG.
    """
    aspects = []
    
    # וודא שיש לפחות שני כוכבים כדי לחשב היבטים
    if len(planet_positions_raw) < 2:
        return aspects

    # נמיר את מילון המיקומים הגולמי לרשימת tuples (שם, קו אורך)
    planets_list = [(name, lon) for name, lon in planet_positions_raw.items()]

    # עובר על כל זוגות הכוכבים האפשריים
    for p1_name, p1_lon in planets_list:
        for p2_name, p2_lon in planets_list:
            # וודא שלא בודקים היבט של כוכב עם עצמו, ושלא בודקים פעמיים (A-B וגם B-A)
            if p1_name >= p2_name: # p1_name >= p2_name for unique pairs
                continue

            # חישוב ההפרש הזוויתי הקצר ביותר
            diff = get_angular_distance(p1_lon, p2_lon)

            # בדיקה מול כל ההיבטים המוגדרים
            for aspect_name, aspect_info in ASPECTS_DEGREE_ORB_CONFIG.items():
                target_degree = aspect_info["degree"]
                orb = aspect_info["orb"]
                aspect_type = aspect_info["type"] # major/soft/minor

                if abs(diff - target_degree) <= orb:
                    # היבט נמצא!
                    actual_orb = round(abs(diff - target_degree), 2)
                    
                    # ניתן להוסיף לוגיקה לחישוב Strength אם נדרש (קרוב יותר למדויק = חזק יותר)
                    strength = "moderate"
                    if actual_orb <= orb * 0.2: # לדוגמה, אם האורב הוא ב-20% מהאורב המקסימלי
                        strength = "strong"
                    elif actual_orb >= orb * 0.8: # אם האורב הוא ב-80% ומעלה מהאורב המקסימלי
                        strength = "weak"

                    aspects.append({
                        "planet1": p1_name,
                        "planet2": p2_name,
                        "aspect": aspect_name,
                        "type": aspect_type,
                        "orb": actual_orb,
                        "exact_degree": target_degree,
                        "strength": strength
                    })
    return aspects


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
                "planets": sorted([p1_name, p2_name, p3_name]), # מיון כדי להבטיח עקביות
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
        return None # לא מספיק כוכבים ליצירת מגן דוד

    all_trines = find_all_grand_trines(positions)
    
    # חפש זוגות של Grand Trines
    for gt1, gt2 in itertools.combinations(all_trines, 2):
        # וודא שהכוכבים בשני הטרינים שונים לחלוטין
        planets_gt1 = set(gt1["planets"])
        planets_gt2 = set(gt2["planets"])
        
        if len(planets_gt1.union(planets_gt2)) != 6:
            continue # יש כוכבים חופפים או יותר מ-6 כוכבים ייחודיים

        # יש לנו 6 כוכבים ייחודיים המחולקים לשני טרינים.
        # כעת נבדוק אם הם יוצרים 3 אופוזיציות
        sod_planets = list(planets_gt1.union(planets_gt2))
        
        # ניצור רשימה של כל הפלנטות והמיקומים שלהן מתוך ה-6 שנבחרו
        current_sod_positions = {p_name: positions[p_name] for p_name in sod_planets}
        current_sod_planets_list = [(p, current_sod_positions[p]) for p in sod_planets]

        found_oppositions = []
        # נבדוק אופוזיציות בין כל זוג אפשרי בתוך קבוצת 6 הכוכבים
        for p_a, p_b in itertools.combinations(current_sod_planets_list, 2):
            dist = get_angular_distance(p_a[1], p_b[1])
            if abs(dist - 180) <= GEOMETRIC_PATTERN_ORB:
                found_oppositions.append(tuple(sorted([p_a[0], p_b[0]]))) # מיון שמות כדי למנוע כפילויות

        # מגן דוד דורש בדיוק 3 אופוזיציות ייחודיות (כל כוכב משתתף באופוזיציה אחת)
        # דרך פשוטה לוודא: מספר האופוזיציות הוא 3, ואין כוכבים שחוזרים על עצמם באופוזיציות
        unique_opposition_planets = set()
        for op in found_oppositions:
            unique_opposition_planets.add(op[0])
            unique_opposition_planets.add(op[1])
        
        if len(found_oppositions) == 3 and len(unique_opposition_planets) == 6:
            # וודא ששתי קבוצות הכוכבים (של הטרינים) משתלבות ליצירת אופוזיציות
            # מורכב יותר, אך ניתן לבדוק על ידי כך שכל פלנטה בטרין אחד מוצאת את בן זוגה באופוזיציה בטרין השני
            # לצורך מימוש ראשוני, נסתפק בכך שנמצאו 3 אופוזיציות נפרדות בין 6 הפלנטות הללו.

            return {
                "name": "Star of David",
                "planets": sod_planets,
                "trines": [gt1, gt2],
                "oppositions": found_oppositions
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


# --- חדש: פונקציית עזר לזיהוי בית עבור כוכב נתון ---
def get_house_for_planet(planet_longitude: float, house_cusps_longitudes: List[float]) -> int:
    """
    מזהה באיזה בית נמצא כוכב לכת מסוים, בהתבסס על קו האורך שלו
    ועל קווי האורך של נקודות הבתים (cusps).
    הפונקציה מקבלת רשימה ממוינת של קווי אורך של הבתים.
    """
    # וודא שיש 12 קווי אורך של בתים
    if len(house_cusps_longitudes) < 12:
        # print(f"Warning: Not enough house cusps provided ({len(house_cusps_longitudes)} instead of 12).") # הסר הערה אם רוצים לוגינג
        return -1 # ערך שגיאה או לא ידוע

    # pyswisseph.houses() מחזיר 12 קווי אורך של בתים (cusps[0]...cusps[11]).
    # בית 1 מתחיל מ-cusps[0], בית 2 מ-cusps[1] וכו'.
    # הכוכב נמצא בבית X אם קו האורך שלו בין cusp[X-1] ל-cusp[X].
    # יש לטפל במקרה של מעבר דרך 360/0 מעלות (נקודת האפס במזל טלה).

    for i in range(12):
        house_start_lon = house_cusps_longitudes[i]
        house_end_lon = house_cusps_longitudes[(i + 1) % 12] # הבית הבא, חזרה ל-0 אחרי בית 12

        # טיפול במעבר דרך נקודת האפס (0/360 מעלות)
        if house_start_lon <= house_end_lon:
            # מקרה רגיל: בית בתוך טווח רציף
            if house_start_lon <= planet_longitude < house_end_lon:
                return i + 1
        else:
            # מקרה שבו הבית חוצה את נקודת האפס (לדוגמה, בית 12 חוצה 0 מעלות טלה)
            if planet_longitude >= house_start_lon or planet_longitude < house_end_lon:
                return i + 1
    
    return -1 # אם לא נמצא בית (שזה לא אמור לקרות אם כל הבתים מכוסים)


def calculate_birth_chart(year: int, month: int, day: int, lat: float, long_geo: float) -> Dict[str, Any]:
    """
    מחשב מפת לידה מלאה: מיקומי כוכבים, מיקומי בתים, תבניות גאומטריות והיבטים.
    מקבל שנה בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    מחזירה מילון מובנה הכולל את כל הנתונים.
    **מעודכן**: כולל מהירויות כוכבים, תבניות גאומטריות והיבטים בפלט.
    **חדש**: כולל מספר בית לכל כוכב.
    """
    # se.utc_to_jd מקבל שנה כפי שהיא, כולל שליליות עבור שנים אסטרונומיות.
    # שנה 0 היא שנת 1 לפנה"ס. שנה -1 היא שנת 2 לפנה"ס וכן הלאה.
    jd_utc = se.utc_to_jd(year, month, day, 0, 0, 0)[1]

    # יצירת אובייקט תאריך עבור שמירה/הצגה (ייתכן ותהיה פחות עמידה לשנים מאוד שליליות ב-Python מובנה)
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
        "all_aspects": [], # חדש: רשימה של כל ההיבטים
        "detected_patterns": [] 
    }

    raw_planet_longitudes = {} # נשמור קווי אורך כדי להעביר לפונקציות זיהוי תבניות והיבטים
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
    
    # חישוב מיקומי בתים וקווי האורך הגולמיים של הבתים (cusps)
    # se.houses מחזיר 12 קווי אורך עבור הבתים ו-2 עבור ה-ASC/MC
    raw_cusps, ascmc = se.houses(jd_utc, lat, long_geo, ord('P')) # Placidus
    raw_cusps_longitudes = [raw_cusps[i] for i in range(12)] # קווי אורך של 12 הבתים

    calculated_house_cusps_dict = _calculate_houses(jd_utc, lat, long_geo) # משתמש ב-raw_cusps_longitudes
    chart_data["house_positions"] = calculated_house_cusps_dict

    # --- חדש: זיהוי בית לכל כוכב לכת ---
    for planet_name_hebrew, details in chart_data["planet_positions"].items():
        if "error" not in details:
            planet_lon = details["longitude"]
            house_number = get_house_for_planet(planet_lon, raw_cusps_longitudes)
            details["house"] = house_number # הוספת מספר הבית לכוכב

    # --- חדש: חישוב כל ההיבטים ---
    chart_data["all_aspects"] = calculate_all_aspects(raw_planet_longitudes)


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
    
    # בדיקת היבטים
    print("\n--- בדיקת היבטים בחישוב מפת לידה לדוגמה (1986-04-13) ---")
    if "all_aspects" in chart and chart["all_aspects"]:
        for aspect in chart["all_aspects"]:
            print(f"  - {aspect['planet1']}-{aspect['planet2']}: {aspect['aspect']} ({aspect['orb']} מעלות, סוג {aspect['type']}) - {aspect['strength']}")
    else:
        print("  - לא נמצאו היבטים.")

    # בדיקת מספרי בתים
    print("\n--- בדיקת מספרי בתים לכוכבים (1986-04-13) ---")
    if "planet_positions" in chart:
        for planet, details in chart["planet_positions"].items():
            if "house" in details:
                print(f"  - {planet}: בית {details['house']}")


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
