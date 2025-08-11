import swisseph as se
from datetime import date, timedelta, datetime
from typing import List, Dict, Any, Optional
import itertools

# מילון שמות כוכבים
PLANETS_HEBREW_MAP = {
    "שמש": se.SUN, "ירח": se.MOON, "מרקורי": se.MERCURY, "ונוס": se.VENUS,
    "מאדים": se.MARS, "יופיטר": se.JUPITER, "שבתאי": se.SATURN,
    "אורנוס": se.URANUS, "נפטון": se.NEPTUNE, "פלוטו": se.PLUTO,
    "ראש תלי": se.TRUE_NODE, "זנב תלי": se.MEAN_NODE # ניתן להשתמש ב-MEAN_NODE או ב-TRUE_NODE
}

# מילון היבטים עם מעלות ואורבים
ASPECTS_HEBREW_MAP = {
    "צמידות": {"degree": 0, "orb": 2},
    "היפוך": {"degree": 180, "orb": 2},
    "טרין": {"degree": 120, "orb": 2},
    "סקסטיל": {"degree": 60, "orb": 2},
    "ריבוע": {"degree": 90, "orb": 2}
}

# הגדרת אורב כללי לתבניות גאומטריות
GEOMETRIC_PATTERN_ORB = 5 # מעלות

def calculate_next_date_for_graph(current_date: date, resolution: str) -> date:
    """מחשבת את התאריך הבא בהתאם לרזולוציה עבור גרף סינוסי.
    מקבלת אובייקט date.
    """
    if resolution == "weekly":
        return current_date + timedelta(days=7)
    elif resolution == "monthly":
        # לטיפול בחודשים עם מספר ימים שונה
        current_month = current_date.month
        current_year = current_date.year
        current_day = current_date.day

        # ננסה לקדם בחודש. אם התאריך לא חוקי, נקדם ליום האחרון של החודש הבא.
        try:
            next_month = current_month + 1
            next_year = current_year
            if next_month > 12:
                next_month = 1
                next_year += 1
            
            # בדיקה עבור שנים אסטרונומיות
            if next_year == 0: # אם השנה היא 0 בפורמט אסטרונומי, זהו 1 לפנה"ס
                # אובייקט date לא תומך בשנה 0, אז נטפל בזה ידנית או נמנע מיצירה ישירה
                # עבור תאריך הבא, נדלג לשנה הבאה אם אנחנו צפויים להגיע ל-0.
                # במקרה של חישוב תאריך הבא, סביר ש-current_date כבר תקין
                pass 
            
            next_date_candidate = date(next_year, next_month, current_day)
            return next_date_candidate
        except ValueError:
            # אם היום לא חוקי בחודש הבא (לדוגמה, 31 בפברואר), קח את היום האחרון של החודש הבא
            next_month = current_month + 1
            next_year = current_year
            if next_month > 12:
                next_month = 1
                next_year += 1

            # יצירת היום הראשון של החודש הבא, ואז חיסור יום אחד
            first_day_of_next_month = date(next_year, next_month, 1)
            return first_day_of_next_month - timedelta(days=1)

    elif resolution == "yearly":
        return date(current_date.year + 1, current_date.month, current_date.day)
    else: # daily או כל דבר אחר כברירת מחדל
        return current_date + timedelta(days=1)

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

        # בדיקה האם כל ההפרשים קרובים ל-120 מעלות
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


def generate_sine_chart_data(
    planet_names: List[str], 
    start_year: int, start_month: int, start_day: int,
    end_year: int, end_month: int, end_day: int,
    aspects_to_find: Optional[List[str]] = None, 
    resolution: str = "daily", 
    patterns_to_find: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """
    מייצר נתונים אסטרולוגיים עבור גרף סינוסי, כולל מיקומי כוכבים, היבטים ותבניות גאומטריות,
    לפי רזולוציה.
    מקבל שנות התחלה וסיום בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    """
    se.set_ephe_path('ephe')
    
    graph_data = []
    
    # יצירת אובייקטי date מתאריכי התחלה וסיום
    # אובייקט date ב-Python לא תומך בשנת 0000, אז נשתמש ב-datetime ונטפל בהמרות
    # Swiss Ephemeris (se.utc_to_jd) כן תומך בשנה 0.
    # כאן אנו מייצרים אובייקט date עבור הלוגיקה של timedelta,
    # כאשר שנת 0 מומרת ל-1 לספירה, ועל כן יש לקחת בחשבון את ההבדל בהצגה
    # מול ה-frontend שידע להציג נכון את ה-BCE.
    start_dt = datetime(start_year, start_month, start_day)
    end_dt = datetime(end_year, end_month, end_day)

    current_date = start_dt.date() # שימוש באובייקט date
    
    planet_ids = {name: PLANETS_HEBREW_MAP[name] for name in planet_names if name in PLANETS_HEBREW_MAP}
    
    # חיפוש היבטים
    target_aspect_degrees = {ASPECTS_HEBREW_MAP[aspect]["degree"]: ASPECTS_HEBREW_MAP[aspect]["orb"] for aspect in aspects_to_find or [] if aspect in ASPECTS_HEBREW_MAP}

    print(f"מתחיל יצירת נתונים לגרף סינוסי מ- {start_dt.strftime('%Y-%m-%d')} עד {end_dt.strftime('%Y-%m-%d')} ברזולוציית {resolution}...")

    while current_date <= end_dt.date():
        # se.utc_to_jd מקבל שנה, חודש, יום בנפרד (כולל שנים אסטרונומיות 0 או שליליות)
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        daily_positions = {}
        for planet_name, planet_id in planet_ids.items():
            try:
                planet_pos, _ = se.calc_ut(jd_utc, planet_id)
                daily_positions[planet_name] = round(planet_pos[0], 2)
            except Exception as e:
                # הדפסת שגיאה אם חישוב כוכב נכשל
                print(f"שגיאה בחישוב מיקום {planet_name} בתאריך {current_date}: {e}")
                daily_positions[planet_name] = None # סמן כוכב שלא חושב כראוי

        daily_aspects = []
        # אם יש יותר מכוכב אחד ברשימה וכדאי לבדוק היבטים
        if len(planet_names) >= 2 and target_aspect_degrees:
            # נבדוק רק בין שני הכוכבים הראשונים ברשימה שנבחרו לגרף, או שנקבל פרמטרים ספציפיים
            # NOTE: if we want to check aspects for ALL pairs of selected planets, this needs to be a nested loop
            # for p1_name, p2_name in itertools.combinations(planet_names, 2): # אם רוצים לכל הזוגות
            #     if (p1_name in daily_positions and daily_positions[p1_name] is not None and
            #         p2_name in daily_positions and daily_positions[p2_name] is not None):
            #         lon1 = daily_positions[p1_name]
            #         lon2 = daily_positions[p2_name]
            #         ... (שאר לוגיקת ההיבט)
            # מכיוון שה-frontend מאפשר רק 2 כוכבים לתרשים גלי רגיל (אלא אם נבחרו תבניות),
            # נתמקד בבדיקת היבטים בין שני הכוכבים הראשונים.
            if (planet_names[0] in daily_positions and daily_positions[planet_names[0]] is not None and
                planet_names[1] in daily_positions and daily_positions[planet_names[1]] is not None):

                lon1 = daily_positions[planet_names[0]]
                lon2 = daily_positions[planet_names[1]]
                
                diff = abs(lon1 - lon2)
                if diff > 180:
                    diff = 360 - diff

                for aspect_degree, orb_val in target_aspect_degrees.items():
                    if abs(diff - aspect_degree) <= orb_val:
                        # מציאת שם ההיבט מהמילון המקורי
                        aspect_name_found = next((name for name, info in ASPECTS_HEBREW_MAP.items() if info["degree"] == aspect_degree), "היבט לא ידוע")
                        daily_aspects.append({
                            "planet1": planet_names[0],
                            "planet2": planet_names[1],
                            "aspect": aspect_name_found,
                            "orb": round(abs(diff - aspect_degree), 2)
                        })
        
        # זיהוי תבניות גאומטריות
        daily_patterns = []
        if patterns_to_find:
            if "Grand Trine" in patterns_to_find:
                if len(daily_positions) >= 3:
                    found_gts = find_all_grand_trines(daily_positions)
                    if found_gts:
                        daily_patterns.extend(found_gts)
            
            if "Star of David" in patterns_to_find:
                if len(daily_positions) >= 6: # מגן דוד דורש 6 כוכבים לפחות
                    sod = find_star_of_david(daily_positions)
                    if sod:
                        daily_patterns.append(sod)


        graph_data.append({
            "date": str(current_date), # נשמור את התאריך בפורמט YYYY-MM-DD
            "positions": daily_positions,
            "aspects": daily_aspects,
            "patterns": daily_patterns
        })
        
        current_date = calculate_next_date_for_graph(current_date, resolution)
    
    print("יצירת נתונים לגרף סינוסי הסתיימה.")
    return graph_data

if __name__ == "__main__":
    # בדיקת הפונקציה
    
    # בדיקה עבור שנים חיוביות
    print("--- יצירת נתונים לגרף סינוסי עבור שמש, ירח, מאדים (עם Grand Trine) (2023) ---")
    sine_data_gt = generate_sine_chart_data(
        ["שמש", "ירח", "מאדים"], 
        2023, 1, 1, 2023, 1, 31, 
        ["צמידות", "היפוך"], 
        resolution="daily", 
        patterns_to_find=["Grand Trine"]
    )
    print(f"נמצאו {len(sine_data_gt)} נקודות נתונים.")
    for item in sine_data_gt[:5]:
        print(f" - {item}")
    if len(sine_data_gt) > 5:
        print("...")

    print("\n--- יצירת נתונים לגרף סינוסי עבור יופיטר, שבתאי (ללא Grand Trine) (2023) ---")
    sine_data_no_pattern = generate_sine_chart_data(
        ["יופיטר", "שבתאי"], 
        2023, 1, 1, 2023, 1, 31, 
        resolution="weekly"
    )
    print(f"נמצאו {len(sine_data_no_pattern)} נקודות נתונים.")
    for item in sine_data_no_pattern[:5]:
        print(f" - {item}")
    if len(sine_data_no_pattern) > 5:
        print("...")

    print("\n--- יצירת נתונים לגרף סינוסי עם בדיקת מגן דוד (טווח רחב יותר, כולל שנים לפנה''ס) ---")
    # 1980-01-01 עד 2000-01-01 - שנים חיוביות
    start_sod_test_positive = 1980
    end_sod_test_positive = 2000
    sod_planets_test = ["שמש", "ירח", "מאדים", "יופיטר", "שבתאי", "נפטון"]
    
    sine_data_sod_positive = generate_sine_chart_data(
        sod_planets_test, 
        start_sod_test_positive, 1, 1, 
        end_sod_test_positive, 1, 1, 
        resolution="monthly", 
        patterns_to_find=["Star of David"]
    )
    
    print(f"נמצאו {len(sine_data_sod_positive)} נקודות נתונים.")
    sod_found_count_positive = 0
    for item in sine_data_sod_positive:
        if item["patterns"]:
            print(f" - {item['date']}: {item['patterns']}")
            sod_found_count_positive += 1
    print(f"נמצאו {sod_found_count_positive} תאריכים עם תבניות מגן דוד (שנים חיוביות).")


    # בדיקה עבור שנים לפני הספירה
    print("\n--- יצירת נתונים לגרף סינוסי עם בדיקת מגן דוד (טווח לפני הספירה: 1563 לפנה''ס - 1500 לפנה''ס) ---")
    # 1563 לפנה"ס היא שנה -1562 בפורמט אסטרונומי
    # 1500 לפנה"ס היא שנה -1499 בפורמט אסטרונומי
    start_sod_test_bce_year = -1562 # 1563 BCE
    end_sod_test_bce_year = -1499 # 1500 BCE

    sine_data_sod_bce = generate_sine_chart_data(
        sod_planets_test, 
        start_sod_test_bce_year, 1, 1, 
        end_sod_test_bce_year, 1, 1, 
        resolution="monthly", 
        patterns_to_find=["Star of David"]
    )
    
    print(f"נמצאו {len(sine_data_sod_bce)} נקודות נתונים.")
    sod_found_count_bce = 0
    for item in sine_data_sod_bce:
        if item["patterns"]:
            print(f" - {item['date']}: {item['patterns']}")
            sod_found_count_bce += 1
    print(f"נמצאו {sod_found_count_bce} תאריכים עם תבניות מגן דוד (שנים לפני הספירה).")

