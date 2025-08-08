import json
import os
from datetime import datetime

LIBRARY_FILE = "library.json"

def save_search(tool_name, title, search_params, results, tags=None):
    """
    שומר חיפוש חדש בקובץ הספרייה.
    
    פרמטרים:
    tool_name (str): שם הכלי שבו בוצע החיפוש.
    title (str): כותרת לחיפוש.
    search_params (dict): פרמטרים של החיפוש.
    results (list): תוצאות החיפוש.
    tags (list, optional): רשימת תגיות.
    """
    if tags is None:
        tags = []
        
    try:
        with open(LIBRARY_FILE, 'r', encoding='utf-8') as f:
            library = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        library = []

    new_search = {
        "id": len(library) + 1,
        "tool": tool_name,
        "title": title,
        "search_params": search_params,
        "results": results,
        "tags": tags,
        "timestamp": datetime.now().isoformat()
    }
    
    library.append(new_search)
    
    with open(LIBRARY_FILE, 'w', encoding='utf-8') as f:
        json.dump(library, f, indent=2, ensure_ascii=False)
    
    print(f"חיפוש נשמר בהצלחה בספרייה. ID: {new_search['id']}")

def get_all_searches():
    """מחזיר את כל החיפושים השמורים בספרייה."""
    try:
        with open(LIBRARY_FILE, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return []

def search_library(query=None, by_tool=None, by_tag=None):
    """
    מחפש בספרייה לפי מלל חופשי, כלי או תגית.
    """
    all_searches = get_all_searches()
    if not all_searches:
        return []

    results = []
    for search in all_searches:
        match = True
        
        # חיפוש לפי כלי
        if by_tool and search["tool"] != by_tool:
            match = False
        
        # חיפוש לפי תגית
        if by_tag and by_tag not in search["tags"]:
            match = False
            
        # חיפוש מלל חופשי (בכותרת או בפרמטרים)
        if query and query.lower() not in search["title"].lower() and query.lower() not in json.dumps(search["search_params"]).lower():
            match = False
        
        if match:
            results.append(search)
            
    return results

def delete_search(search_id):
    """מוחק חיפוש מהספרייה על פי ה-ID שלו."""
    all_searches = get_all_searches()
    updated_library = [s for s in all_searches if s["id"] != search_id]
    
    with open(LIBRARY_FILE, 'w', encoding='utf-8') as f:
        json.dump(updated_library, f, indent=2, ensure_ascii=False)
        
    print(f"חיפוש עם ID {search_id} נמחק בהצלחה.")

# --- קוד לבדיקת הפונקציה ---
if __name__ == "__main__":
    print("--- בדיקת מודול הספרייה ---")

    # בדיקה 1: שמירת חיפוש לדוגמה
    sample_search_data = {
        "tool": "historical_pattern_finder",
        "title": "צמידות יופיטר-שבתאי",
        "search_params": {"planet1": "יופיטר", "planet2": "שבתאי", "aspect": "צמידות"},
        "results": [{"date": "2020-12-21", "orb": 0.1}],
        "tags": ["יופיטר", "שבתאי", "היבטים"]
    }
    save_search(
        sample_search_data["tool"],
        sample_search_data["title"],
        sample_search_data["search_params"],
        sample_search_data["results"],
        sample_search_data["tags"]
    )
    
    # בדיקה 2: חיפוש בספרייה
    print("\nכל החיפושים השמורים:")
    all_searches = get_all_searches()
    print(json.dumps(all_searches, indent=2, ensure_ascii=False))

    # בדיקה 3: חיפוש לפי תגית
    print("\nחיפוש לפי תגית 'יופיטר':")
    jupiter_searches = search_library(by_tag="יופיטר")
    print(json.dumps(jupiter_searches, indent=2, ensure_ascii=False))

    # בדיקה 4: מחיקת חיפוש
    if all_searches:
        search_id_to_delete = all_searches[0]["id"]
        delete_search(search_id_to_delete)
        print("\nספרייה לאחר מחיקה:")
        print(json.dumps(get_all_searches(), indent=2, ensure_ascii=False))