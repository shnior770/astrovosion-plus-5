import sys
from datetime import date
from prettytable import PrettyTable
import json
from historical_pattern_finder import find_historical_pattern, find_historical_aspect, PLANETS_HEBREW_MAP, SIGNS_HEBREW_MAP, ASPECTS_HEBREW_MAP
from sine_graph_generator import generate_sine_graph_data
from library_manager import save_search, get_all_searches, search_library, delete_search
from astro import calculate_birth_chart

def show_birth_chart(target_date=None, lat=None, lon=None, data_to_display=None):
    """
    מבקש תאריך, קו רוחב וקו אורך ומציג את מפת הלידה (כוכבים ובתים) בטבלה.
    """
    if target_date is None:
        print("\n" + "="*50)
        print("        הצגת גלגל המזלות לתאריך נבחר")
        print("="*50)
        try:
            year = int(input("הכנס שנה (לדוגמה, 1990): "))
            month = int(input("הכנס חודש (1-12): "))
            day = int(input("הכנס יום (1-31): "))
            lat = float(input("הכנס קו רוחב (לדוגמה, 31.7683): "))
            lon = float(input("הכנס קו אורך (לדוגמה, 35.2137): "))
            target_date = date(year, month, day)
        except ValueError:
            print("קלט לא חוקי. אנא ודא שהכנסת מספרים תקינים.")
            return

    if data_to_display is None:
        print(f"\nמתחיל חישוב מפת לידה עבור {target_date}...")
        full_chart = calculate_birth_chart(target_date, lat, lon)
        data_to_display = full_chart

    if data_to_display:
        # הצגת כוכבים
        table_planets = PrettyTable()
        table_planets.field_names = ["כוכב", "מזל", "מעלה"]
        for planet, data in data_to_display.items():
            if planet != "בתים":
                table_planets.add_row([planet, data["sign"], round(data["degree_in_sign"], 2)])
        
        print("\n" + "="*50)
        print(f"תוצאות עבור {target_date}:")
        print("="*50)
        print("מיקומי כוכבים:")
        print(table_planets)

        # הצגת בתים
        table_houses = PrettyTable()
        table_houses.field_names = ["בית", "מזל", "מעלה"]
        for house, data in data_to_display["בתים"].items():
            table_houses.add_row([house, data["sign"], round(data["degree_in_sign"], 2)])
        
        print("\nמיקומי בתים:")
        print(table_houses)

        save_option = input("האם ברצונך לשמור חיפוש זה בספרייה? (כן/לא): ")
        if save_option.lower() == 'כן':
            title = input("הכנס כותרת לחיפוש זה: ")
            save_search("גלגל מזלות", title, {"date": str(target_date), "lat": lat, "lon": lon}, data_to_display)
    else:
        print("לא ניתן לחשב או להציג את המיקומים לתאריך זה.")


def find_historical_events(search_params=None, data_to_display=None):
    """
    מאפשר למשתמש לבחור בין חיפוש תבנית (כוכב במזל) או חיפוש היבט.
    """
    if search_params is None:
        print("\n" + "="*50)
        print("        חיפוש אירועים היסטוריים")
        print("בחר אפשרות:")
        print("1. חיפוש תבנית (כוכב במזל ובמעלה)")
        print("2. חיפוש היבט (בין שני כוכבים)")
        print("3. חזרה לתפריט הראשי")
        print("="*50)
        try:
            choice = int(input("הכנס את מספר הבחירה: "))
        except ValueError:
            print("קלט לא חוקי. אנא הכנס מספר בין 1 ל-3.")
            return

        if choice == 1:
            print("\n--- חיפוש תבנית ---")
            print(f"כוכבים זמינים: {list(PLANETS_HEBREW_MAP.keys())}")
            planet_name = input("הכנס שם כוכב: ")
            print(f"מזלות זמינים: {list(SIGNS_HEBREW_MAP.values())}")
            sign_name = input("הכנס שם מזל: ")
            try:
                start_year = int(input("הכנס שנת התחלה: "))
                end_year = int(input("הכנס שנת סיום: "))
                search_params = {"type": "pattern", "planet": planet_name, "sign": sign_name, "start_year": start_year, "end_year": end_year}
            except ValueError:
                print("קלט לא חוקי.")
                return
        elif choice == 2:
            print("\n--- חיפוש היבט ---")
            print(f"כוכבים זמינים: {list(PLANETS_HEBREW_MAP.keys())}")
            planet1 = input("הכנס שם כוכב ראשון: ")
            planet2 = input("הכנס שם כוכב שני: ")
            print(f"היבטים זמינים: {list(ASPECTS_HEBREW_MAP.keys())}")
            aspect_name = input("הכנס שם היבט: ")
            try:
                start_year = int(input("הכנס שנת התחלה: "))
                end_year = int(input("הכנס שנת סיום: "))
                search_params = {"type": "aspect", "planet1": planet1, "planet2": planet2, "aspect": aspect_name, "start_year": start_year, "end_year": end_year}
            except ValueError:
                print("קלט לא חוקי.")
                return
        elif choice == 3:
            return
        else:
            print("בחירה לא חוקית. אנא נסה שוב.")
            return
    
    if data_to_display is None:
        start_date = date(search_params["start_year"], 1, 1)
        end_date = date(search_params["end_year"], 12, 31)
        
        if search_params["type"] == "pattern":
            results = find_historical_pattern(search_params["planet"], search_params["sign"], start_date, end_date)
            data_to_display = results
            table_fields = ["תאריך", "כוכב", "מזל", "מעלה"]
            table_rows = [[row["date"], search_params["planet"], row["sign"], round(row["degree_in_sign"], 2)] for row in results[:10]]
        elif search_params["type"] == "aspect":
            results = find_historical_aspect(search_params["planet1"], search_params["planet2"], search_params["aspect"], start_date, end_date)
            data_to_display = results
            table_fields = ["תאריך", "כוכבים", "היבט", "אורב"]
            table_rows = [[row["date"], f"{row['planet1']}-{row['planet2']}", row["aspect"], round(row["orb"], 2)] for row in results[:10]]
    else:
        if search_params["type"] == "pattern":
            table_fields = ["תאריך", "כוכב", "מזל", "מעלה"]
            table_rows = [[row["date"], search_params["planet"], row["sign"], round(row["degree_in_sign"], 2)] for row in data_to_display[:10]]
        else:
            table_fields = ["תאריך", "כוכבים", "היבט", "אורב"]
            table_rows = [[row["date"], f"{row['planet1']}-{row['planet2']}", row["aspect"], round(row["orb"], 2)] for row in data_to_display[:10]]


    if data_to_display:
        table = PrettyTable()
        table.field_names = table_fields
        for row in table_rows:
            table.add_row(row)
        print(table)

        if search_params is not None and not data_to_display:
            save_option = input("האם ברצונך לשמור חיפוש זה בספרייה? (כן/לא): ")
            if save_option.lower() == 'כן':
                title = input("הכנס כותרת לחיפוש זה: ")
                tags_input = input("הכנס תגיות מופרדות בפסיק (לדוגמה: יופיטר,צמידות): ")
                tags = [t.strip() for t in tags_input.split(',')]
                save_search("חיפוש היסטורי", title, search_params, data_to_display, tags)
    else:
        print("לא נמצאו תוצאות בטווח התאריכים.")

def show_sine_graph(search_params=None, data_to_display=None):
    """
    מבקש מהמשתמש לבחור כוכבים, טווח תאריכים והיבטים, ומציג נתונים גרפיים.
    """
    if search_params is None:
        print("\n" + "="*50)
        print("        הצגת גרף סינוסי")
        print("="*50)
        
        print(f"כוכבים זמינים: {list(PLANETS_HEBREW_MAP.keys())}")
        planets_input = input("הכנס שמות כוכבים מופרדים בפסיק (לדוגמה: שמש,ירח): ")
        planet_names = [p.strip() for p in planets_input.split(',')]
        
        print(f"היבטים זמינים: {list(ASPECTS_HEBREW_MAP.keys())}")
        aspects_input = input("הכנס שמות היבטים לחיפוש (מופרדים בפסיק, לדוגמה: צמידות,מולות): ")
        aspects_to_find = [a.strip() for a in aspects_input.split(',')] if aspects_input else None
        
        try:
            start_year = int(input("הכנס שנת התחלה: "))
            end_year = int(input("הכנס שנת סיום: "))
            search_params = {"planets": planet_names, "aspects": aspects_to_find, "start_year": start_year, "end_year": end_year}
        except ValueError:
            print("קלט לא חוקי. אנא הכנס מספרים תקינים.")
            return

    if data_to_display is None:
        start_date = date(search_params["start_year"], 1, 1)
        end_date = date(search_params["end_year"], 12, 31)
        graph_data = generate_sine_graph_data(search_params["planets"], start_date, end_date, search_params["aspects"])
        data_to_display = graph_data
    
    if data_to_display:
        print("\n" + "="*50)
        print("נתונים גרפיים (מוצגים בפורמט JSON):")
        print(json.dumps(data_to_display, indent=2, ensure_ascii=False))

        if search_params is not None:
            save_option = input("האם ברצונך לשמור חיפוש זה בספרייה? (כן/לא): ")
            if save_option.lower() == 'כן':
                title = input("הכנס כותרת לחיפוש זה: ")
                tags_input = input("הכנס תגיות מופרדות בפסיק (לדוגמה: שמש,גרף): ")
                tags = [t.strip() for t in tags_input.split(',')]
                save_search("גרף סינוסי", title, search_params, data_to_display, tags)
    else:
        print("לא ניתן לחשב או להציג את נתוני הגרף.")

def manage_library():
    """
    מנהל את הספרייה ומאפשר צפייה, חיפוש, מחיקה וטעינת חיפושים שמורים.
    """
    while True:
        print("\n" + "="*50)
        print("        ניהול ספרייה")
        print("בחר אפשרות:")
        print("1. הצג את כל החיפושים השמורים")
        print("2. חפש בספרייה")
        print("3. מחק חיפוש מהספרייה")
        print("4. טען חיפוש מהספרייה")
        print("5. חזור לתפריט הראשי")
        print("="*50)

        try:
            choice = int(input("הכנס את מספר הבחירה: "))
        except ValueError:
            print("קלט לא חוקי. אנא הכנס מספר בין 1 ל-5.")
            continue
        
        if choice == 1:
            all_searches = get_all_searches()
            if not all_searches:
                print("הספרייה ריקה.")
            else:
                table = PrettyTable()
                table.field_names = ["ID", "כותרת", "כלי", "תגיות", "תאריך שמירה"]
                for row in all_searches:
                    table.add_row([row["id"], row["title"], row["tool"], ", ".join(row["tags"]), row["timestamp"][:10]])
                print(table)
        
        elif choice == 2:
            query = input("הכנס מילת מפתח לחיפוש (השאר ריק אם אין): ")
            by_tool = input("הכנס שם כלי לחיפוש (השאר ריק אם אין): ")
            by_tag = input("הכנס תגית לחיפוש (השאר ריק אם אין): ")
            
            results = search_library(query if query else None, by_tool if by_tool else None, by_tag if by_tag else None)
            if not results:
                print("לא נמצאו תוצאות תואמות.")
            else:
                table = PrettyTable()
                table.field_names = ["ID", "כותרת", "כלי", "תגיות", "תאריך שמירה"]
                for row in results:
                    table.add_row([row["id"], row["title"], row["tool"], ", ".join(row["tags"]), row["timestamp"][:10]])
                print(table)
        
        elif choice == 3:
            try:
                search_id = int(input("הכנס את ה-ID של החיפוש למחיקה: "))
                delete_search(search_id)
            except ValueError:
                print("קלט לא חוקי. ה-ID חייב להיות מספר.")
        
        elif choice == 4:
            try:
                search_id = int(input("הכנס את ה-ID של החיפוש לטעינה: "))
                load_and_display_search(search_id)
            except ValueError:
                print("קלט לא חוקי. ה-ID חייב להיות מספר.")

        elif choice == 5:
            break
        else:
            print("בחירה לא חוקית. אנא נסה שוב.")

def load_and_display_search(search_id):
    """
    טוען חיפוש מהספרייה ומציג אותו מחדש.
    """
    all_searches = get_all_searches()
    search = next((s for s in all_searches if s["id"] == search_id), None)
    
    if search is None:
        print(f"לא נמצא חיפוש עם ID: {search_id}")
        return
        
    print(f"\nטוען את החיפוש '{search['title']}'...")
    
    if search["tool"] == "גלגל מזלות":
        target_date = date.fromisoformat(search["search_params"]["date"])
        lat = search["search_params"]["lat"]
        lon = search["search_params"]["lon"]
        show_birth_chart(target_date=target_date, lat=lat, lon=lon, data_to_display=search["results"])
    elif search["tool"] == "חיפוש היסטורי":
        find_historical_events(search_params=search["search_params"], data_to_display=search["results"])
    elif search["tool"] == "גרף סינוסי":
        show_sine_graph(search_params=search["search_params"], data_to_display=search["results"])
    else:
        print(f"הכלי '{search['tool']}' אינו נתמך לטעינה.")


def main_menu():
    """
    מציג את התפריט הראשי ומקבל בחירה מהמשתמש.
    """
    print("\n" + "="*50)
    print("ברוכים הבאים למערכת 'אסטרוויז'ן'")
    print("בחר אפשרות מהרשימה:")
    print("1. הצג גלגל מזלות")
    print("2. חיפוש אירועים היסטוריים")
    print("3. הצג גרף סינוסי")
    print("4. ניהול ספרייה")
    print("5. יציאה")
    print("="*50)

    try:
        choice = int(input("הכנס את מספר הבחירה: "))
        return choice
    except ValueError:
        print("קלט לא חוקי. אנא הכנס מספר בין 1 ל-5.")
        return None

def main():
    while True:
        choice = main_menu()
        if choice == 1:
            show_birth_chart()
        elif choice == 2:
            find_historical_events()
        elif choice == 3:
            show_sine_graph()
        elif choice == 4:
            manage_library()
        elif choice == 5:
            print("היציאה מהמערכת...")
            sys.exit()
        else:
            print("בחירה לא חוקית. אנא נסה שוב.")

if __name__ == "__main__":
    main()