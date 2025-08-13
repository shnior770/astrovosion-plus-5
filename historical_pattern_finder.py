import pyswisseph as se
from datetime import date, timedelta, datetime
from typing import List, Dict, Any, Optional

# ייבוא נתונים ופונקציות מודול הליבה astro.py
# זה מבטיח עקביות ומניעת כפילויות בהגדרות
from astro import (
    PLANETS_HEBREW_MAP,
    SIGNS_HEBREW_MAP_DEGREE_RANGES, # נשתמש במילון המזלות עם טווחי המעלות
    get_sign_from_longitude,
    GEOMETRIC_PATTERN_ORB, # נייבא את האורב הגאומטרי הקבוע
    get_angular_distance,
    find_all_grand_trines,
    find_star_of_david,
    calculate_birth_chart # ייבוא calculate_birth_chart לצורך 'all_planets_at_date'
)

# הגדרת נתיב לקבצי האפמריס (חובה עבור swisseph)
se.set_ephe_path('ephe')

# תקופות מסלוליות ממוצעות של כוכבים (בימים, לצורך הערכה)
# המידע הזה ישמש לאופטימיזציית חיפושים ארוכים
PLANETARY_ORBITAL_PERIODS_DAYS = {
    "שמש": 365.25,
    "ירח": 27.32,
    "מרקורי": 87.97,
    "ונוס": 224.7,
    "מאדים": 686.98,
    "יופיטר": 4332.59,   # ~11.86 שנים
    "שבתאי": 10759.22,  # ~29.46 שנים
    "אורנוס": 30688.5,   # ~84.02 שנים
    "נפטון": 60182.0,   # ~164.79 שנים
    "פלוטו": 90560.0,   # ~247.92 שנים
    "ראש תלי": 6798.3,   # ~18.6 שנים (מחזור נסיגה)
    "זנב תלי": 6798.3    # זהה לראש תלי
}

# הגדרת סף לזיהוי כוכב כ"איטי" (לדוגמה, תקופה מעל 10 שנים)
SLOW_PLANET_THRESHOLD_DAYS = 365 * 10 

def calculate_next_date_for_search(current_date: date, resolution: str) -> date:
    """מחשבת את התאריך הבא בהתאם לרזולוציה.
    מקבלת אובייקט date.
    """
    if resolution == "weekly":
        return current_date + timedelta(days=7)
    elif resolution == "monthly":
        year = current_date.year
        month = current_date.month + 1
        if month > 12:
            month = 1
            year += 1
        # מנסה לשמור על היום בחודש, או לעבור לסוף החודש אם הוא לא קיים (למשל, 31 בפברואר)
        try:
            # שימוש ב-datetime לצורך טיפול בשנה 0 במידת הצורך
            return date(year, month, current_date.day)
        except ValueError:
            # אם היום לא חוקי בחודש הבא (לדוגמה, 31 בפברואר), קח את היום האחרון של החודש הבא
            next_month = current_date.month + 1
            next_year = current_date.year
            if next_month > 12:
                next_month = 1
                next_year += 1

            # יצירת היום הראשון של החודש הבא, ואז חיסור יום אחד
            first_day_of_next_month = date(next_year, next_month, 1)
            return first_day_of_next_month - timedelta(days=1)
    elif resolution == "yearly":
        return date(current_date.year + 1, current_date.month, current_date.day)
    else: # daily
        return current_date + timedelta(days=1)

def find_constellation_planet_in_sign(
    planet_name: str, 
    sign_name: str, 
    start_year: int, start_month: int, start_day: int,
    end_year: int, end_month: int, end_day: int,
    resolution: str = "daily", 
    degree: Optional[int] = None,
    max_results: Optional[int] = None # חדש: פרמטר להגבלת תוצאות
) -> List[Dict[str, Any]]:
    """
    מחפש תאריכים שבהם כוכב מסוים נמצא במזל ספציפי, ואולי גם במעלה ספציפית.
    מקבל שנות התחלה וסיום בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    **מעודכן**: תומך ב-max_results וכולל את כל מיקומי הכוכבים בתאריך שנמצא.
    """
    se.set_ephe_path('ephe')
    
    results = []
    
    # יצירת אובייקטי date מתאריכי התחלה וסיום
    start_dt = datetime(start_year, start_month, start_day)
    end_dt = datetime(end_year, end_month, end_day)

    current_date = start_dt.date() # שימוש באובייקט date

    planet_id = PLANETS_HEBREW_MAP.get(planet_name)
    if planet_id is None:
        # החזרת שגיאה בפורמט שצד הלקוח יוכל לזהות כבעיה בפרמטר קלט
        return [{"error_type": "InvalidInput", "message": f"שם כוכב לא חוקי: {planet_name}"}]

    # ודא שהמזל קיים במילון טווחי המעלות
    if sign_name not in SIGNS_HEBREW_MAP_DEGREE_RANGES:
        return [{"error_type": "InvalidInput", "message": f"שם מזל לא חוקי: {sign_name}"}]
    
    sign_start_deg, sign_end_deg = SIGNS_HEBREW_MAP_DEGREE_RANGES[sign_name]

    print(f"מתחיל חיפוש היסטורי עבור {planet_name} ב{sign_name} (מעלה {degree if degree is not None else 'כלשהי'}) מ- {start_dt.strftime('%Y-%m-%d')} עד {end_dt.strftime('%Y-%m-%d')} ברזולוציית {resolution}...")

    while current_date <= end_dt.date():
        # עצירה אם הגענו למגבלת התוצאות
        if max_results is not None and len(results) >= max_results:
            break

        # se.utc_to_jd מקבל שנה, חודש, יום בנפרד (כולל שנים אסטרונומיות 0 או שליליות)
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        try:
            planet_pos, _ = se.calc_ut(jd_utc, planet_id)
            longitude = planet_pos[0] % 360 # ודא שהמעלות בין 0 ל-360

            # בדיקה אם הכוכב נמצא במזל הנכון
            if sign_start_deg <= longitude < sign_end_deg:
                degree_in_sign = int(longitude) % 30
                # אם נתון גם מעלה, בדוק התאמה
                if degree is None or degree_in_sign == degree:
                    # חדש: קריאה ל-calculate_birth_chart כדי לקבל את כל מיקומי הכוכבים
                    # נצטרך לספק ל-calculate_birth_chart את השנה בפורמט הקלנדרי המקורי
                    # אבל calculate_birth_chart כבר מטפלת בהמרת השנה האסטרונומית
                    # ולכן נשתמש בשנה האסטרונומית ישירות עבור calculate_birth_chart
                    
                    # זיהוי השנה האסטרונומית של current_date
                    astro_year = current_date.year
                    if current_date.year <= 0: # אם זה 1 לפנה"ס או לפני
                        astro_year = -(current_date.year - 1) # המרה לשנה קלנדרית רגילה
                        if astro_year == 1: # 1 לפנה"ס
                            astro_year = 0 # חזרה לפורמט אסטרונומי
                        elif astro_year > 1:
                            astro_year = - (astro_year -1)

                    
                    # קריאה ל-calculate_birth_chart כדי לקבל את כל מיקומי הכוכבים ליום זה
                    # נצטרך לספק קו רוחב וקו אורך דמה כי find_constellation_planet_in_sign לא מקבל אותם.
                    # נניח ש-lat, long_geo הם פרמטרים נוספים אם נרצה חישוב מדויק של כל המפה,
                    # אך לצורך 'all_planets_at_date' המפתח ביקש רק מיקומים
                    
                    # לטובת מילוי 'all_planets_at_date', נבצע חישוב מיקומי כוכבים בלבד
                    # אם נצטרך את כל המידע כמו ב-calculate_birth_chart, נצטרך להעביר את lat, long
                    
                    # פשוט נחשב את מיקומי כל הכוכבים ליום הזה
                    all_planets_daily_pos = {}
                    for p_name_hebrew, p_id in PLANETS_HEBREW_MAP.items():
                        try:
                            daily_planet_pos, _ = se.calc_ut(jd_utc, p_id, se.SWE_FLG_SPEED)
                            all_planets_daily_pos[p_name_hebrew] = {
                                "longitude": round(daily_planet_pos[0] % 360, 2),
                                "speed": round(daily_planet_pos[3], 4),
                                "is_retrograde": daily_planet_pos[3] < 0
                            }
                        except Exception as inner_e:
                            print(f"שגיאה בחישוב כל הכוכבים עבור {current_date}: {inner_e}")
                            all_planets_daily_pos[p_name_hebrew] = {"error": "Failed to calculate position"}


                    results.append({
                        "date": str(current_date),
                        "planet_name": planet_name,
                        "sign_name": sign_name,
                        "longitude": round(longitude, 2),
                        "degree_in_sign": degree_in_sign,
                        "all_planets_at_date": all_planets_daily_pos # חדש: כל מיקומי הכוכבים
                    })
        except Exception as e:
            print(f"שגיאה בחישוב מיקום {planet_name} בתאריך {current_date}: {e}")
            # נמשיך הלאה גם אם יש שגיאה ביום ספציפי

        current_date = calculate_next_date_for_search(current_date, resolution)
    
    print("חיפוש היסטורי הסתיים.")
    return results

def find_constellation_aspect(
    planet1_name: str, planet2_name: str, aspect_name: str, 
    start_year: int, start_month: int, start_day: int,
    end_year: int, end_month: int, end_day: int,
    resolution: str = "daily", orb: float = 1.0, # אורב ברירת מחדל
    max_results: Optional[int] = None # חדש: פרמטר להגבלת תוצאות
) -> List[Dict[str, Any]]:
    """
    מחפש תאריכים שבהם נוצר היבט ספציפי בין שני כוכבים.
    מקבל שנות התחלה וסיום בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).
    **מעודכן**: תומך ב-max_results וכולל את כל מיקומי הכוכבים בתאריך שנמצא.
    """
    se.set_ephe_path('ephe')
    
    results = []
    
    # יצירת אובייקטי date מתאריכי התחלה וסיום
    start_dt = datetime(start_year, start_month, start_day)
    end_dt = datetime(end_year, end_month, end_day)

    current_date = start_dt.date() # שימוש באובייקט date

    planet1_id = PLANETS_HEBREW_MAP.get(planet1_name)
    planet2_id = PLANETS_HEBREW_MAP.get(planet2_name)
    
    if planet1_id is None or planet2_id is None:
        return [{"error_type": "InvalidInput", "message": f"שם כוכב לא חוקי: {planet1_name} או {planet2_name}"}]

    # השתמש במילון היבטים המפורט מ-astro.py
    from astro import ASPECTS_DEGREE_ORB_CONFIG # ייבוא חדש!
    
    aspect_info = ASPECTS_DEGREE_ORB_CONFIG.get(aspect_name)
    if aspect_info is None:
        return [{"error_type": "InvalidInput", "message": f"שם היבט לא חוקי: {aspect_name}"}]
    
    target_degree = aspect_info["degree"]
    # השתמש באורב המוגדר מראש במילון ASPECT_DEGREE_ORB_CONFIG, אלא אם הועבר אורב ספציפי
    # נגדיר את זה כברירת מחדל ואפשרות ל-override
    effective_orb = orb if orb is not None else aspect_info["orb"] 


    print(f"מתחיל חיפוש היסטורי עבור היבט {aspect_name} בין {planet1_name} ל-{planet2_name} מ- {start_dt.strftime('%Y-%m-%d')} עד {end_dt.strftime('%Y-%m-%d')} ברזולוציית {resolution}...")

    while current_date <= end_dt.date():
        # עצירה אם הגענו למגבלת התוצאות
        if max_results is not None and len(results) >= max_results:
            break

        # se.utc_to_jd מקבל שנה, חודש, יום בנפרד (כולל שנים אסטרונומיות 0 או שליליות)
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        try:
            pos1, _ = se.calc_ut(jd_utc, planet1_id)
            pos2, _ = se.calc_ut(jd_utc, planet2_id)
            
            lon1 = pos1[0] % 360
            lon2 = pos2[0] % 360

            diff = abs(lon1 - lon2)
            if diff > 180:
                diff = 360 - diff

            if abs(diff - target_degree) <= effective_orb: # שימוש ב-effective_orb
                # חדש: קריאה לחישוב כל הכוכבים ליום זה
                all_planets_daily_pos = {}
                for p_name_hebrew, p_id in PLANETS_HEBREW_MAP.items():
                    try:
                        daily_planet_pos, _ = se.calc_ut(jd_utc, p_id, se.SWE_FLG_SPEED)
                        all_planets_daily_pos[p_name_hebrew] = {
                            "longitude": round(daily_planet_pos[0] % 360, 2),
                            "speed": round(daily_planet_pos[3], 4),
                            "is_retrograde": daily_planet_pos[3] < 0
                        }
                    except Exception as inner_e:
                        print(f"שגיאה בחישוב כל הכוכבים עבור {current_date}: {inner_e}")
                        all_planets_daily_pos[p_name_hebrew] = {"error": "Failed to calculate position"}

                results.append({
                    "date": str(current_date),
                    "planet1": planet1_name,
                    "planet2": planet2_name,
                    "aspect": aspect_name,
                    "actual_difference": round(diff, 2),
                    "orb": round(abs(diff - target_degree), 2),
                    "all_planets_at_date": all_planets_daily_pos # חדש: כל מיקומי הכוכבים
                })
        except Exception as e:
            print(f"שגיאה בחישוב מיקום בתאריך {current_date}: {e}")

        current_date = calculate_next_date_for_search(current_date, resolution)
    
    print("חיפוש היסטורי היבטים הסתיים.")
    return results

def find_complex_constellation(
    conditions: List[Dict[str, Any]], 
    start_year: int, start_month: int, start_day: int,
    end_year: int, end_month: int, end_day: int,
    resolution: str = "daily",
    max_results: Optional[int] = None # חדש: פרמטר להגבלת תוצאות
) -> List[Dict[str, Any]]:
    """
    מחפש תאריכים שבהם קבוצת כוכבים עומדת במספר תנאים בו-זמנית.
    לדוגמה: שמש בטלה 15 מעלות, וירח בסרטן.
    מקבל שנות התחלה וסיום בפורמט אסטרונומי (שלילי עבור שנים לפני הספירה).

    בנוסף, הפונקציה כוללת אופטימיזציית חישוב:
    אם החיפוש מערב כוכבים איטיים (כמו פלוטו) ברזולוציה יומית/שבועית על פני טווח ארוך,
    תתקבל אזהרה לגבי יעילות.
    **מעודכן**: תומך ב-max_results וכולל את כל מיקומי הכוכבים בתאריך שנמצא.
    """
    se.set_ephe_path('ephe')
    
    results = []
    
    # יצירת אובייקטי date מתאריכי התחלה וסיום
    start_dt = datetime(start_year, start_month, start_day)
    end_dt = datetime(end_year, end_month, end_day)

    current_date = start_dt.date() # שימוש באובייקט date

    # --- אופטימיזציית חישוב: זיהוי כוכבים איטיים ---
    involved_planets = [cond['planet_name'] for cond in conditions]
    slowest_period_days = 0
    are_slow_planets_involved = False
    for planet_name in involved_planets:
        period = PLANETARY_ORBITAL_PERIODS_DAYS.get(planet_name, 0) # 0 אם לא נמצא (לדוגמה: בתים)
        if period > SLOW_PLANET_THRESHOLD_DAYS:
            are_slow_planets_involved = True
        if period > slowest_period_days:
            slowest_period_days = period

    # הערכת משך החיפוש הכולל בימים
    total_search_days = (end_dt.date() - start_dt.date()).days # השתמש ב-date objects להפרש

    # אזהרה על חיפוש לא יעיל בכוכבים איטיים
    if are_slow_planets_involved:
        if total_search_days > 365 * 2 and (resolution == "daily" or resolution == "weekly"):
            optimization_warning = {
                "type": "warning",
                "message": (
                    f"אזהרה: חיפוש תבנית זו מערב כוכבים איטיים (כמו {involved_planets[0] if are_slow_planets_involved else ''}). "
                    f"חיפוש ברזולוציית '{resolution}' לאורך {total_search_days} ימים עלול להיות לא יעיל. "
                    f"שקול רזולוציה 'חודשית' או 'שנתית' לחיפוש מהיר יותר."
                )
            }
            results.append(optimization_warning) # הוסף את האזהרה לתוצאות (ניתן גם לזרוק חריגה או לשלוח בלוג)

    print(f"מתחיל חיפוש מורכב עבור {len(conditions)} תנאים מ- {start_dt.strftime('%Y-%m-%d')} עד {end_dt.strftime('%Y-%m-%d')} ברזולוציית {resolution}...")

    # כעת, נחזור על הלוגיקה המקורית לחיפוש תבניות מורכבות:
    while current_date <= end_dt.date():
        # עצירה אם הגענו למגבלת התוצאות
        if max_results is not None and len(results) >= max_results:
            break

        # se.utc_to_jd מקבל שנה, חודש, יום בנפרד (כולל שנים אסטרונומיות 0 או שליליות)
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        
        all_conditions_met = True
        daily_positions_for_conditions = {} # נאסוף את המיקומים עבור התנאים הספציפיים
        
        for condition in conditions:
            planet_name = condition.get('planet_name')
            sign_name = condition.get('sign_name')
            target_degree = condition.get('degree')

            planet_id = PLANETS_HEBREW_MAP.get(planet_name)
            if planet_id is None:
                all_conditions_met = False
                break # כוכב לא חוקי, דילוג על היום
            
            # ודא שהמזל קיים במילון טווחי המעלות
            if sign_name not in SIGNS_HEBREW_MAP_DEGREE_RANGES:
                all_conditions_met = False
                break # מזל לא חוקי, דילוג על היום

            sign_start_deg, sign_end_deg = SIGNS_HEBREW_MAP_DEGREE_RANGES[sign_name]

            try:
                planet_pos, _ = se.calc_ut(jd_utc, planet_id)
                longitude = planet_pos[0] % 360

                if not (sign_start_deg <= longitude < sign_end_deg):
                    all_conditions_met = False
                    break # התנאי לכוכב זה לא מתקיים
                
                degree_in_sign = int(longitude) % 30
                if target_degree is not None and degree_in_sign != target_degree:
                    all_conditions_met = False
                    break # התנאי למעלה לא מתקיים

                daily_positions_for_conditions[planet_name] = {
                    "longitude": round(longitude, 2),
                    "sign": sign_name,
                    "degree_in_sign": degree_in_sign
                }

            except Exception as e:
                print(f"שגיאה בחישוב מיקום {planet_name} בתאריך {current_date}: {e}")
                all_conditions_met = False
                break # שגיאה, דילוג על היום

        if all_conditions_met:
            # חדש: קריאה לחישוב כל הכוכבים ליום זה
            all_planets_daily_pos = {}
            for p_name_hebrew, p_id in PLANETS_HEBREW_MAP.items():
                try:
                    daily_planet_pos, _ = se.calc_ut(jd_utc, p_id, se.SWE_FLG_SPEED)
                    all_planets_daily_pos[p_name_hebrew] = {
                        "longitude": round(daily_planet_pos[0] % 360, 2),
                        "speed": round(daily_planet_pos[3], 4),
                        "is_retrograde": daily_planet_pos[3] < 0
                    }
                except Exception as inner_e:
                    print(f"שגיאה בחישוב כל הכוכבים עבור {current_date}: {inner_e}")
                    all_planets_daily_pos[p_name_hebrew] = {"error": "Failed to calculate position"}

            results.append({
                "date": str(current_date),
                "planets_meeting_conditions": daily_positions_for_conditions,
                "all_planets_at_date": all_planets_daily_pos # חדש: כל מיקומי הכוכבים
            })

        current_date = calculate_next_date_for_search(current_date, resolution)
    
    print("חיפוש מורכב הסתיים.")
    return results


if __name__ == "__main__":
    # בדיקת פונקציות (מופעלות רק כשהקובץ מורץ ישירות)
    
    # בדיקת find_constellation_planet_in_sign
    print("--- בדיקת find_constellation_planet_in_sign (שמש בטלה 15, 2023) ---")
    results1 = find_constellation_planet_in_sign("שמש", "טלה", 2023, 3, 1, 2023, 4, 30, "daily", 15)
    print(f"נמצאו {len(results1)} תאריכים עבור שמש בטלה 15: {results1}")

    # בדיקה עם max_results
    print("\n--- בדיקת find_constellation_planet_in_sign (עם max_results=2) ---")
    results_max = find_constellation_planet_in_sign("שמש", "טלה", 2023, 3, 1, 2023, 4, 30, "daily", max_results=2)
    print(f"נמצאו {len(results_max)} תאריכים עבור שמש בטלה (max_results=2): {results_max}")


    print("\n--- בדיקת find_constellation_planet_in_sign (שמש בטלה 15, 1563 לפנה''ס - 1562 לפנה''ס) ---")
    # 1563 לפנה"ס היא שנה -1562 בפורמט אסטרונומי
    # 1562 לפנה"ס היא שנה -1561 בפורמט אסטרונומי
    results_bce_planet = find_constellation_planet_in_sign("שמש", "טלה", -1562, 3, 1, -1561, 4, 30, "daily", 15)
    print(f"נמצאו {len(results_bce_planet)} תאריכים עבור שמש בטלה 15 (לפנה''ס): {results_bce_planet}")


    # בדיקת find_constellation_aspect
    print("\n--- בדיקת find_constellation_aspect (צמידות ירח-שמש, 2023) ---")
    results2 = find_constellation_aspect("ירח", "שמש", "צמידות", 2023, 1, 1, 2023, 3, 1, "daily")
    print(f"נמצאו {len(results2)} תאריכים עבור צמידות ירח-שמש: {results2}")

    # בדיקה עם max_results
    print("\n--- בדיקת find_constellation_aspect (עם max_results=2) ---")
    results_max_aspect = find_constellation_aspect("ירח", "שמש", "צמידות", 2023, 1, 1, 2023, 3, 1, "daily", max_results=2)
    print(f"נמצאו {len(results_max_aspect)} תאריכים עבור צמידות ירח-שמש (max_results=2): {results_max_aspect}")


    print("\n--- בדיקת find_constellation_aspect (צמידות ירח-שמש, 1 לפנה''ס - 1 לספירה) ---")
    # 1 לפנה"ס היא שנה 0 בפורמט אסטרונומי
    # 1 לספירה היא שנה 1 בפורמט אסטרונומי
    results_bce_aspect = find_constellation_aspect("ירח", "שמש", "צמידות", 0, 1, 1, 1, 1, 1, "daily")
    print(f"נמצאו {len(results_bce_aspect)} תאריכים עבור צמידות ירח-שמש (לפנה''ס): {results_bce_aspect}")


    # בדיקת find_complex_constellation עם אופטימיזציה (כוכב איטי)
    print("\n--- בדיקת find_complex_constellation (יופיטר במזל ספציפי, טווח ארוך, 2000-2005) ---")
    complex_cond_slow = [
        {"planet_name": "יופיטר", "sign_name": "שור", "degree": 1} # יופיטר כוכב איטי
    ]
    results_complex_slow = find_complex_constellation(complex_cond_slow, 2000, 1, 1, 2005, 1, 1, "daily")
    print(f"נמצאו {len(results_complex_slow)} תאריכים עבור חיפוש מורכב (יופיטר):")
    for r in results_complex_slow:
        print(r)

    # בדיקה עם max_results
    print("\n--- בדיקת find_complex_constellation (עם max_results=2) ---")
    results_max_complex = find_complex_constellation(complex_cond_slow, 2000, 1, 1, 2005, 1, 1, "daily", max_results=2)
    print(f"נמצאו {len(results_max_complex)} תאריכים עבור חיפוש מורכב (max_results=2):")
    for r in results_max_complex:
        print(r)


    print("\n--- בדיקת find_complex_constellation (שמש וירח במזל ספציפי, טווח קצר, 2024) ---")
    complex_cond_fast = [
        {"planet_name": "שמש", "sign_name": "דלי", "degree": 5},
        {"planet_name": "ירח", "sign_name": "טלה", "degree": 10}
    ]
    results_complex_fast = find_complex_constellation(complex_cond_fast, 2024, 1, 1, 2024, 2, 28, "daily")
    print(f"נמצאו {len(results_complex_fast)} תאריכים עבור חיפוש מורכב (שמש וירח):")
    for r in results_complex_fast:
        print(r)

    print("\n--- בדיקת find_complex_constellation (שמש וירח במזל ספציפי, 1 לפנה''ס) ---")
    # 1 לפנה"ס היא שנה 0 אסטרונומית
    complex_cond_bce = [
        {"planet_name": "שמש", "sign_name": "טלה", "degree": 15},
        {"planet_name": "ירח", "sign_name": "מאזניים", "degree": 15}
    ]
    results_complex_bce = find_complex_constellation(complex_cond_bce, 0, 1, 1, 0, 12, 31, "daily")
    print(f"נמצאו {len(results_complex_bce)} תאריכים עבור חיפוש מורכב (שמש וירח, לפנה''ס):")
    for r in results_complex_bce:
        print(r)
