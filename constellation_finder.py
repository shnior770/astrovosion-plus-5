import swisseph as se
from datetime import date, timedelta
import json

# מילון שמות כוכבים
PLANETS_HEBREW_MAP = {
    "שמש": se.SUN, "ירח": se.MOON, "מרקורי": se.MERCURY, "ונוס": se.VENUS,
    "מאדים": se.MARS, "יופיטר": se.JUPITER, "שבתאי": se.SATURN,
    "אורנוס": se.URANUS, "נפטון": se.NEPTUNE, "פלוטו": se.PLUTO
}

# מילון קונסטלציות אסטרונומיות וטווח המעלות שלהן
# (הערה: אלה טווחים משוערים ויכולים להשתנות, ניתן לשפר אותם בעתיד)
CONSTELLATIONS_MAP = {
    "טלה": (0, 30),
    "שור": (30, 60),
    "תאומים": (60, 90),
    "סרטן": (90, 120),
    "אריה": (120, 150),
    "בתולה": (150, 180),
    "מאזניים": (180, 210),
    "עקרב": (210, 240),
    "קשת": (240, 270),
    "גדי": (270, 300),
    "דלי": (300, 330),
    "דגים": (330, 360)
}

def get_constellation_name(longitude):
    """
    מחזירה את שם הקונסטלציה על פי קו האורך.
    """
    for name, (start_lon, end_lon) in CONSTELLATIONS_MAP.items():
        if start_lon <= longitude < end_lon:
            return name
    return "קונסטלציה לא ידועה"

def find_constellation(planet_name, start_date, end_date):
    """
    סורקת תקופה היסטורית ומחזירה את הקונסטלציה של כוכב מסוים בכל יום.
    """
    planet_id = PLANETS_HEBREW_MAP.get(planet_name)
    if planet_id is None:
        print(f"שם כוכב לא חוקי: {planet_name}")
        return []

    results = []
    current_date = start_date
    delta = timedelta(days=1)
    
    print(f"מתחיל סריקה של הקונסטלציה עבור {planet_name} מ- {start_date} עד {end_date}...")
    
    while current_date <= end_date:
        jd_utc = se.utc_to_jd(current_date.year, current_date.month, current_date.day, 0, 0, 0)[1]
        planet_pos, _ = se.calc_ut(jd_utc, planet_id)
        longitude = planet_pos[0]
        
        constellation = get_constellation_name(longitude)
        
        results.append({
            "date": str(current_date),
            "constellation": constellation,
            "longitude": round(longitude, 2)
        })
        
        current_date += delta
        
    print("סריקה הסתיימה.")
    return results

# --- קוד לבדיקת הפונקציה ---
if __name__ == "__main__":
    # הגדרת טווח קצר לבדיקות מהירות
    start_date = date(2025, 1, 1)
    end_date = date(2025, 1, 31)
    
    # לולאה על כל כוכבי הלכת במילון
    for planet_name in PLANETS_HEBREW_MAP.keys():
        print("\n" + "="*50 + "\n")
        print(f"--- חיפוש: מיקום הקונסטלציה עבור {planet_name} ---")
        
        constellations = find_constellation(planet_name, start_date, end_date)
        
        if constellations:
            print(f"נמצאו {len(constellations)} נקודות נתונים.")
            # מדפיס רק את 5 התוצאות הראשונות כדי למנוע הדפסה ארוכה
            print(json.dumps(constellations[:5], indent=2, ensure_ascii=False))
            
            # אם יש יותר מ-5 תוצאות, ציין זאת
            if len(constellations) > 5:
                print("...")
        else:
            print(f"לא נמצאו אירועים בטווח התאריכים הקצר עבור {planet_name}.")