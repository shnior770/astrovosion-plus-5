import pyswisseph as se
from datetime import date, timedelta, datetime
from typing import List, Dict, Any, Optional
import itertools

# הגדרת נתיב לקבצי האפמריס (חובה עבור swisseph)
se.set_ephe_path('ephe')

# מילון שמות כוכבים
PLANETS_HEBREW_MAP = {
    "שמש": se.SUN, "ירח": se.MOON, "מרקורי": se.MERCURY, "ונוס": se.VENUS,
    "מאדים": se.MARS, "יופיטר": se.JUPITER, "שבתאי": se.SATURN,
    "אורנוס": se.URANUS, "נפטון": se.NEPTUNE, "פלוטו": se.PLUTO,
    "ראש תלי": se.TRUE_NODE, "זנב תלי": se.MEAN_NODE # ניתן להשתמש ב-MEAN_NODE או ב-TRUE_NODE
}

# --- עדכון: מילון היבטים עם מעלות ואורבים (כמו ב-astro.py) ---
# אנו מייבאים את זה מ-astro.py כדי להבטיח עקביות בכל המערכת
from astro import ASPECTS_DEGREE_ORB_CONFIG, GEOMETRIC_PATTERN_ORB
# שימו לב: ASPECTS_HEBREW_MAP המקורי הוחלף ב-ASPECTS_DEGREE_ORB_CONFIG מ-astro.py
# ו-GEOMETRIC_PATTERN_ORB מיובא גם כן.


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

        try:
            next_month = current_month + 1
            next_year = current_year
            if next_month > 12:
                next_month = 1
                next_year += 1
            
            # עבור שנים אסטרונומיות (שליליות או 0), אובייקט date עשוי להיות בעייתי.
            # אנו מניחים ש-datetime.date(year, month, day) יטפל בזה באופן סביר
            # עבור הטווחים הנפוצים, ושהשנה כבר הומרה לפורמט אסטרונומי ב-app.py.
            return date(next_year, next_month, current_day)
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

# --- חדש: פונקציה לזיהוי הצטלבות מדויקת בין שני כוכבים (לצורך גרף סינוסי) ---
def find_exact_intersections_and_aspects_between_points(
    p1_name: str, p1_lon1: float, p1_lon2: float,
    p2_name: str, p2_lon1: float, p2_lon2: float,
    current_date: date, next_date: date,
    aspects_to_find_config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """
    מבצע אינטרפולציה לינארית בין שתי נקודות זמן (יומיות/נקודות רזולוציה)
    כדי לזהות הצטלבויות מדויקות או היבטים מדויקים (אורב 0) בין שני כוכבים.
    """
    found_events = []

    # חישוב שינוי מעלות ליום (או לטווח הזמן הזה)
    delta_p1 = p1_lon2 - p1_lon1
    delta_p2 = p2_lon2 - p2_lon1
    
    # אם כוכב חוצה 360/0, נתאים את הדלתא
    if abs(delta_p1) > 180: delta_p1 -= 360 * (delta_p1 / abs(delta_p1))
    if abs(delta_p2) > 180: delta_p2 -= 360 * (delta_p2 / abs(delta_p2))

    # ההפרש היחסי במעלה בין שני הכוכבים
    relative_speed = delta_p1 - delta_p2

    # אם אין תנועה יחסית משמעותית, אין הצטלבות או היבט מדויק
    if abs(relative_speed) < 0.0001: # סף קטן מאוד למניעת חלוקה באפס
        return found_events

    # מספר הימים/צעדים בין current_date ל-next_date
    num_steps = (next_date - current_date).days # לרוב יהיה 1 ליומי, 7 לשבועי וכו'

    # בדיקת הצטלבות (Conjunction) - כאשר ההפרש בין הכוכבים הוא 0 מעלות
    initial_diff = get_angular_distance(p1_lon1, p2_lon1)
    
    # חיפוש היבטים מדויקים
    for aspect_name, aspect_info in ASPECTS_DEGREE_ORB_CONFIG.items():
        target_degree = aspect_info["degree"]
        
        # חישוב ההפרש בין ההיבט הקיים ליום הנוכחי לבין היעד
        # נחשב את ההפרש היחסי (distance to target aspect)
        current_aspect_diff_from_target = get_angular_distance(p1_lon1, p2_lon1) - target_degree
        next_aspect_diff_from_target = get_angular_distance(p1_lon2, p2_lon2) - target_degree
        
        # כדי לבדוק חציית אפס, אנו רוצים לראות אם הסימנים של ההפרשים שונים
        # כלומר, אם ההיבט עבר דרך נקודת היעד (target_degree)
        
        # טיפול מיוחד בנקודת המעבר מ-360 ל-0 מעלות (או להיפך)
        # זה קצת מורכב יותר מאינטרפולציה לינארית פשוטה
        # לצורך פשטות ודיוק, נסתכל על הקו הישר בין שתי הנקודות
        
        # הדרך הפשוטה ביותר לזיהוי חצייה:
        # אם יש היבט ביום הראשון (עם אורב קטן, נניח 0.1)
        # ואם יש היבט ביום השני (עם אורב קטן, נניח 0.1)
        # זה מצביע על מעבר קרוב.
        # אינטרפולציה לינארית בסיסית:
        
        # חיפוש חצייה של היבט מדויק (אורב 0)
        # ננסה למצוא את הנקודה בה (lon1 + speed1*t) - (lon2 + speed2*t) = target_degree +/- 360*N
        # (lon1-lon2) + (speed1-speed2)*t = target_degree +/- 360*N
        # (initial_diff) + (relative_speed)*t = target_degree +/- 360*N

        # פתרון עבור t (זמן): t = (target_degree - initial_diff + 360*N) / relative_speed
        # נבדוק עבור N= -1, 0, 1
        
        # יש מספר נקודות חצייה אפשריות לאורך 360 מעלות
        # נבדוק 3 מקרים של יעד: target_degree, target_degree+360, target_degree-360
        target_degrees_to_check = [target_degree, target_degree + 360, target_degree - 360]

        for current_target in target_degrees_to_check:
            # חישוב הזמן היחסי (t) בתוך המרווח (0 עד num_steps)
            if relative_speed != 0: # לוודא שאין חלוקה באפס
                t = (current_target - (p1_lon1 - p2_lon1)) / relative_speed
            else:
                continue # אין תנועה יחסית, אין חצייה לינארית

            # אם t נמצא בטווח (0, num_steps), יש חצייה בתוך המרווח הנוכחי
            if 0 <= t <= num_steps: # כולל נקודות הקצה
                # חישוב התאריך המדויק של ההצטלבות/היבט
                exact_date = current_date + timedelta(days=t)
                
                # חישוב מיקום מדויק בזמן ההצטלבות
                exact_p1_lon = (p1_lon1 + delta_p1 * (t / num_steps)) % 360
                exact_p2_lon = (p2_lon1 + delta_p2 * (t / num_steps)) % 360

                found_events.append({
                    "type": aspect_name, # סוג האירוע: שם ההיבט
                    "date": str(exact_date.date()), # תאריך האירוע
                    "exact_time_fraction": t, # חלק היחסי של היום
                    "planets": [p1_name, p2_name],
                    "longitudes_at_event": {
                        p1_name: round(exact_p1_lon, 2),
                        p2_name: round(exact_p2_lon, 2)
                    }
                })
    return found_events


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
    **מעודכן**: כולל זיהוי נקודות הצטלבות/היבטים מדויקים (intersections).
    """
    se.set_ephe_path('ephe')
    
    graph_data = []
    all_intersections = [] # רשימה חדשה לאחסון הצטלבויות מדויקות
    
    # יצירת אובייקטי date מתאריכי התחלה וסיום
    start_dt_obj = datetime(start_year, start_month, start_day).date()
    end_dt_obj = datetime(end_year, end_month, end_day).date()

    current_date = start_dt_obj 
    
    planet_ids = {name: PLANETS_HEBREW_MAP[name] for name in planet_names if name in PLANETS_HEBREW_MAP}
    
    # חיפוש היבטים
    # נשתמש ב-ASPECTS_DEGREE_ORB_CONFIG המיובא מ-astro.py
    target_aspect_degrees_for_daily_check = {
        ASPECTS_DEGREE_ORB_CONFIG[aspect]["degree"]: ASPECTS_DEGREE_ORB_CONFIG[aspect]["orb"] 
        for aspect in aspects_to_find or [] if aspect in ASPECTS_DEGREE_ORB_CONFIG
    }
    
    # הגדרות היבטים מלאות עבור find_exact_intersections_and_aspects_between_points
    # נכין מילון של aspect_name: aspect_info (degree, orb, type)
    aspects_for_exact_intersections = {
        name: info for name, info in ASPECTS_DEGREE_ORB_CONFIG.items() 
        if name in (aspects_to_find or []) # רק היבטים שהמשתמש ביקש
    }


    print(f"מתחיל יצירת נתונים לגרף סינוסי מ- {start_dt_obj.strftime('%Y-%m-%d')} עד {end_dt_obj.strftime('%Y-%m-%d')} ברזולוציית {resolution}...")

    # נשמור את המיקומים מהאיטרציה הקודמת כדי לבצע אינטרפולציה
    previous_daily_positions = {}
    previous_jd_utc = None

    while current_date <= end_dt_obj:
        # se.utc_to_jd מקבל שנה, חודש, יום בנפרד (כולל שנים אסטרונומיות 0 או שליליות)
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        daily_positions = {}
        for planet_name, planet_id in planet_ids.items():
            try:
                # חשוב לקבל גם מהירות לצורך חישוב מדויק יותר של הצטלבויות (אם כי אינטרפולציה לינארית פשוטה משתמשת רק ב-lon1, lon2)
                planet_pos, _ = se.calc_ut(jd_utc, planet_id, se.SWE_FLG_SPEED) 
                daily_positions[planet_name] = round(planet_pos[0], 2)
            except Exception as e:
                print(f"שגיאה בחישוב מיקום {planet_name} בתאריך {current_date}: {e}")
                daily_positions[planet_name] = None 


        # --- זיהוי הצטלבויות (intersections) בין ימים ---
        # נבצע זאת רק אם זו לא האיטרציה הראשונה ורק אם יש לפחות שני כוכבים נבחרים
        if previous_jd_utc is not None and len(planet_ids) >= 2:
            # נבצע בדיקת הצטלבויות לכל זוג כוכבים שנבחר
            selected_planets_for_intersections = list(planet_ids.keys())
            for p1_name, p2_name in itertools.combinations(selected_planets_for_intersections, 2):
                if (p1_name in previous_daily_positions and previous_daily_positions[p1_name] is not None and
                    p2_name in previous_daily_positions and previous_daily_positions[p2_name] is not None and
                    p1_name in daily_positions and daily_positions[p1_name] is not None and
                    p2_name in daily_positions and daily_positions[p2_name] is not None):
                    
                    found_exact_events = find_exact_intersections_and_aspects_between_points(
                        p1_name, previous_daily_positions[p1_name], daily_positions[p1_name],
                        p2_name, previous_daily_positions[p2_name], daily_positions[p2_name],
                        current_date - timedelta(days=(current_date - previous_date).days), # תאריך התחלה של המרווח
                        current_date, # תאריך סיום של המרווח
                        aspects_for_exact_intersections # העבר את הגדרות ההיבטים
                    )
                    if found_exact_events:
                        all_intersections.extend(found_exact_events)


        daily_aspects = []
        # אם יש יותר מכוכב אחד ברשימה וכדאי לבדוק היבטים (עם אורב)
        if len(planet_names) >= 2 and target_aspect_degrees:
            for p1_name, p2_name in itertools.combinations(planet_names, 2):
                if (p1_name in daily_positions and daily_positions[p1_name] is not None and
                    p2_name in daily_positions and daily_positions[p2_name] is not None):

                    lon1 = daily_positions[p1_name]
                    lon2 = daily_positions[p2_name]
                    
                    diff = abs(lon1 - lon2)
                    if diff > 180:
                        diff = 360 - diff

                    for aspect_degree, orb_val in target_aspect_degrees.items():
                        if abs(diff - aspect_degree) <= orb_val:
                            # מציאת שם ההיבט מהמילון המקורי (ASPECTS_DEGREE_ORB_CONFIG)
                            aspect_name_found = next((name for name, info in ASPECTS_DEGREE_ORB_CONFIG.items() if info["degree"] == aspect_degree), "היבט לא ידוע")
                            daily_aspects.append({
                                "planet1": p1_name,
                                "planet2": p2_name,
                                "aspect": aspect_name_found,
                                "orb": round(abs(diff - aspect_degree), 2)
                            })
        
        # זיהוי תבניות גאומטריות
        daily_patterns = []
        if patterns_to_find:
            # נמיר את daily_positions לפורמט של find_all_grand_trines / find_star_of_david
            # (מילון של שם_כוכב: קו_אורך)
            positions_for_patterns = {name: lon for name, lon in daily_positions.items() if lon is not None}

            if "Grand Trine" in patterns_to_find:
                if len(positions_for_patterns) >= 3:
                    found_gts = find_all_grand_trines(positions_for_patterns)
                    if found_gts:
                        daily_patterns.extend(found_gts)
            
            if "Star of David" in patterns_to_find:
                if len(positions_for_patterns) >= 6: # מגן דוד דורש 6 כוכבים לפחות
                    sod = find_star_of_david(positions_for_patterns)
                    if sod:
                        daily_patterns.append(sod)


        graph_data.append({
            "date": str(current_date), # נשמור את התאריך בפורמט YYYY-MM-DD
            "positions": daily_positions,
            "aspects": daily_aspects,
            "patterns": daily_patterns
        })
        
        # שמירת המיקומים הנוכחיים לאיטרציה הבאה
        previous_daily_positions = daily_positions
        previous_jd_utc = jd_utc
        previous_date = current_date # שמור את התאריך הקודם

        current_date = calculate_next_date_for_graph(current_date, resolution)
    
    print("יצירת נתונים לגרף סינוסי הסתיימה.")
    
    return {
        "data_points": graph_data,
        "intersections": all_intersections # חדש: החזרת כל ההצטלבויות המדויקות
    }

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
    print(f"נמצאו {len(sine_data_gt['data_points'])} נקודות נתונים.")
    for item in sine_data_gt['data_points'][:5]:
        print(f" - {item}")
    if len(sine_data_gt['data_points']) > 5:
        print("...")
    print(f"נמצאו {len(sine_data_gt['intersections'])} הצטלבויות מדויקות.")
    for intersection in sine_data_gt['intersections']:
        print(f"  - Intersection: {intersection}")


    print("\n--- יצירת נתונים לגרף סינוסי עבור יופיטר, שבתאי (ללא Grand Trine) (2023) ---")
    sine_data_no_pattern = generate_sine_chart_data(
        ["יופיטר", "שבתאי"], 
        2023, 1, 1, 2023, 1, 31, 
        resolution="weekly"
    )
    print(f"נמצאו {len(sine_data_no_pattern['data_points'])} נקודות נתונים.")
    for item in sine_data_no_pattern['data_points'][:5]:
        print(f" - {item}")
    if len(sine_data_no_pattern['data_points']) > 5:
        print("...")
    print(f"נמצאו {len(sine_data_no_pattern['intersections'])} הצטלבויות מדויקות.")


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
    
    print(f"נמצאו {len(sine_data_sod_positive['data_points'])} נקודות נתונים.")
    sod_found_count_positive = 0
    for item in sine_data_sod_positive['data_points']:
        if item["patterns"]:
            print(f" - {item['date']}: {item['patterns']}")
            sod_found_count_positive += 1
    print(f"נמצאו {sod_found_count_positive} תאריכים עם תבניות מגן דוד (שנים חיוביות).")
    print(f"נמצאו {len(sine_data_sod_positive['intersections'])} הצטלבויות מדויקות.")


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
    
    print(f"נמצאו {len(sine_data_sod_bce['data_points'])} נקודות נתונים.")
    sod_found_count_bce = 0
    for item in sine_data_sod_bce['data_points']:
        if item["patterns"]:
            print(f" - {item['date']}: {item['patterns']}")
            sod_found_count_bce += 1
    print(f"נמצאו {sod_found_count_bce} תאריכים עם תבניות מגן דוד (שנים לפני הספירה).")
    print(f"נמצאו {len(sine_data_sod_bce['intersections'])} הצטלבויות מדויקות.")

