from flask import Flask, render_template, request, jsonify, Response
from datetime import datetime
import os
import sys
import time
import json
import csv # ייבוא מודול csv
from io import StringIO # לייצוא ל-CSV בזיכרון

# ודא שהסקריפט יוכל למצוא את הקבצים
sys.path.append('.')

# ייבוא הפונקציות המעודכנות מהמודולים שלך
try:
    from astro import calculate_birth_chart, GEOMETRIC_PATTERN_ORB, PLANETS_HEBREW_MAP
    # historical_pattern_finder.py ו-sine_graph_generator.py לא סופקו, לכן נניח שייבואן תקין.
    # אם הם לא קיימים/נחוצים לפונקציונליות הבסיסית, ניתן להסיר.
    from historical_pattern_finder import find_constellation_planet_in_sign, find_constellation_aspect, find_complex_constellation
    from sine_graph_generator import generate_sine_chart_data
except ImportError as e:
    print(f"שגיאה בייבוא המודולים. ודא שכל הקבצים נמצאים בתיקייה הנכונה ובמבנה הנכון.")
    print(f"שגיאה: {e}")
    sys.exit(1)


app = Flask(__name__)

# נקודת קצה לדף הבית (הדף היחיד שלנו)
@app.route('/')
def index():
    """טוען את קובץ ה-HTML הראשי המכיל את כל הטאבים."""
    return render_template('index.html')

# פונקציית עזר להמרת שנה קלנדרית לשנה אסטרונומית (תומך לפני הספירה)
def convert_year_to_astronomical(year, is_bce):
    """
    ממיר שנה קלנדרית (למשל, 1 עבור 1 לספירה, 1 עבור 1 לפנה"ס)
    לשנה אסטרונומית (למשל, 1 עבור 1 לספירה, 0 עבור 1 לפנה"ס, -1 עבור 2 לפנה"ס).
    """
    if is_bce:
        # 1 BCE is year 0, 2 BCE is year -1, etc.
        return -(year - 1)
    return year

# נקודת קצה לחישוב מפת לידה/אירוע (קריאה לפונקציה מתוך astro.py)
@app.route('/calculate', methods=['POST'])
def calculate():
    """
    מקבל נתוני תאריך (שנה, חודש, יום, ודגל is_bce), קו רוחב וקו אורך
    ומחשב מפת לידה/אירוע. מחזיר את התוצאה בפורמט JSON.
    **מעודכן**: הפלט כעת כולל מהירויות כוכבים ותבניות גאומטריות מזוהות.
    """
    data = request.json
    year = data.get('year')
    month = data.get('month')
    day = data.get('day')
    is_bce = data.get('is_bce')
    lat = data.get('lat')
    long_geo = data.get('long')
    
    # **חדש/משופר: בדיקות קלט חזקות יותר**
    if any(val is None for val in [year, month, day, is_bce, lat, long_geo]):
        return jsonify({"error": "נתונים חסרים: ודא ששנה, חודש, יום, לפני הספירה, קו רוחב וקו אורך סופקו."}), 400
    
    # ודא שהנתונים הם מהטיפוסים הנכונים
    try:
        year = int(year)
        month = int(month)
        day = int(day)
        is_bce = bool(is_bce) # ודא שזה בוליאני
        lat = float(lat)
        long_geo = float(long_geo)
    except ValueError:
        return jsonify({"error": "שגיאה: שנה, חודש, יום, קו רוחב או קו אורך אינם מספרים חוקיים."}), 400

    # המרת השנה לפורמט אסטרונומי לפני שליחה ל-astro.py
    astronomical_year = convert_year_to_astronomical(year, is_bce)

    try:
        start_time = time.time() # מדידת זמן התחלה
        result = calculate_birth_chart(astronomical_year, month, day, lat, long_geo)
        end_time = time.time() # מדידת זמן סיום
        duration_ms = round((end_time - start_time) * 1000, 2) # משך במילישניות

        response_data = {'chart': result, 'duration_ms': duration_ms, 'geometric_pattern_orb': GEOMETRIC_PATTERN_ORB}
        return jsonify(response_data)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# נקודת קצה חדשה לייצוא נתוני מפה/חיפוש ל-CSV/JSON
@app.route('/export_chart_data', methods=['POST'])
def export_chart_data():
    """
    מקבל נתוני מפה/חיפוש (כמו לנקודת הקצה /calculate או /find_pattern)
    ומחזיר אותם בפורמט CSV או JSON.
    """
    data = request.json
    export_format = request.args.get('format', 'json').lower() # קבלת פורמט מפראמטר URL

    if not data:
        return jsonify({"error": "יש לספק נתונים לייצוא."}), 400

    if export_format == 'json':
        # מחזיר את הנתונים כ-JSON רגיל
        return jsonify(data)
    elif export_format == 'csv':
        # ננסה לייצא נתוני כוכבים כטבלה CSV
        # נניח שהנתונים ב-'chart' -> 'planet_positions'
        # נדרשת לוגיקה חכמה יותר אם רוצים לייצא סוגי נתונים שונים (חיפושים, תבניות)
        if 'chart' in data and 'planet_positions' in data['chart']:
            output = StringIO()
            writer = csv.writer(output)

            # כותב כותרות: Planet, Longitude, Sign, Degree in Sign, Speed, Is Retrograde
            headers = ['Planet', 'Longitude', 'Sign', 'Degree_in_Sign', 'Speed', 'Is_Retrograde']
            writer.writerow(headers)

            # כותב שורות עבור כל כוכב
            for planet_name, details in data['chart']['planet_positions'].items():
                if 'error' not in details: # דלג על כוכבים עם שגיאות
                    row = [
                        planet_name,
                        details.get('longitude'),
                        details.get('sign'),
                        details.get('degree_in_sign'),
                        details.get('speed'),
                        details.get('is_retrograde')
                    ]
                    writer.writerow(row)
            
            return Response(output.getvalue(), mimetype='text/csv', headers={"Content-disposition": "attachment; filename=chart_data.csv"})
        else:
            return jsonify({"error": "פורמט CSV נתמך כרגע רק עבור נתוני מפה עם מיקומי כוכבים."}), 400
    else:
        return jsonify({"error": "פורמט ייצוא לא נתמך. פורמטים נתמכים: 'json', 'csv'."}), 400


# נקודת קצה לביצוע חיפוש כוכב במזל (קריאה לפונקציה find_constellation_planet_in_sign)
@app.route('/find_pattern', methods=['POST'])
def find_pattern():
    """
    מקבל שם כוכב, מזל ותאריכי התחלה וסיום (שנה, חודש, יום, ודגל is_bce)
    כדי למצוא דפוסים היסטוריים. מחזיר את התוצאות בפורמט JSON, כולל זמן ביצוע.
    """
    data = request.json
    planet_name = data.get('planet_name')
    sign_name = data.get('sign_name')
    degree = data.get('degree')
    
    # קבלת רכיבי תאריך התחלה
    start_year = data.get('start_year')
    start_month = data.get('start_month')
    start_day = data.get('start_day')
    is_start_bce = data.get('is_start_bce')

    # קבלת רכיבי תאריך סיום
    end_year = data.get('end_year')
    end_month = data.get('end_month')
    end_day = data.get('end_day')
    is_end_bce = data.get('is_end_bce')

    resolution = data.get('resolution')

    # **חדש/משופר: בדיקות קלט חזקות יותר**
    if any(val is None for val in [planet_name, sign_name, start_year, start_month, start_day, is_start_bce,
                                   end_year, end_month, end_day, is_end_bce, resolution]):
        return jsonify({"error": "נתונים חסרים: ודא שכל שדות התאריך, הכוכב, המזל והרזולוציה סופקו."}), 400

    try:
        start_year_astro = convert_year_to_astronomical(int(start_year), bool(is_start_bce))
        end_year_astro = convert_year_to_astronomical(int(end_year), bool(is_end_bce))

        start_time = time.time()
        result = find_constellation_planet_in_sign(
            planet_name, sign_name, 
            start_year_astro, int(start_month), int(start_day), 
            end_year_astro, int(end_month), int(end_day), 
            resolution, degree if degree is not None else None # העבר None אם ריק
        )
        end_time = time.time()
        duration_ms = round((end_time - start_time) * 1000, 2)

        response_data = {'results': result, 'duration_ms': duration_ms}
        return jsonify(response_data)
    except ValueError as e:
        return jsonify({'error': f"שגיאה בעיבוד תאריך או מעלה: {e}"}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# נקודת קצה לביצוע חיפוש היבטים היסטוריים (קריאה לפונקציה find_constellation_aspect)
@app.route('/find_aspect', methods=['POST'])
def find_aspect():
    """
    מקבל שני שמות כוכבים, שם היבט ותאריכי התחלה וסיום (שנה, חודש, יום, ודגל is_bce)
    כדי למצוא היבטים היסטוריים. מחזיר את התוצאות בפורמט JSON, כולל זמן ביצוע.
    """
    data = request.json
    planet1_name = data.get('planet1_name')
    planet2_name = data.get('planet2_name')
    aspect_name = data.get('aspect_name')
    
    # קבלת רכיבי תאריך התחלה
    start_year = data.get('start_year')
    start_month = data.get('start_month')
    start_day = data.get('start_day')
    is_start_bce = data.get('is_start_bce')

    # קבלת רכיבי תאריך סיום
    end_year = data.get('end_year')
    end_month = data.get('end_month')
    end_day = data.get('end_day')
    is_end_bce = data.get('is_end_bce')

    resolution = data.get('resolution')

    # **חדש/משופר: בדיקות קלט חזקות יותר**
    if any(val is None for val in [planet1_name, planet2_name, aspect_name, start_year, start_month, start_day, is_start_bce,
                                   end_year, end_month, end_day, is_end_bce, resolution]):
        return jsonify({"error": "נתונים חסרים: ודא שכל שדות התאריך, הכוכבים, ההיבט והרזולוציה סופקו."}), 400

    try:
        start_year_astro = convert_year_to_astronomical(int(start_year), bool(is_start_bce))
        end_year_astro = convert_year_to_astronomical(int(end_year), bool(is_end_bce))
        
        start_time = time.time()
        result = find_constellation_aspect(
            planet1_name, planet2_name, aspect_name,
            start_year_astro, int(start_month), int(start_day),
            end_year_astro, int(end_month), int(end_day),
            resolution
        )
        end_time = time.time()
        duration_ms = round((end_time - start_time) * 1000, 2)

        response_data = {'results': result, 'duration_ms': duration_ms}
        return jsonify(response_data)
    except ValueError as e:
        return jsonify({'error': f"שגיאה בעיבוד תאריך: {e}"}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# נקודת קצה חדשה לביצוע חיפוש קונסטלציה מורכבת (קריאה לפונקציה find_complex_constellation)
@app.route('/find_complex_constellation', methods=['POST'])
def find_complex_constellation_route():
    """
    מקבל רשימת תנאים מורכבים ותאריכי התחלה וסיום (שנה, חודש, יום, ודגל is_bce)
    כדי למצוא תבניות מורכבות. מחזיר את התוצאות בפורמט JSON, כולל זמן ביצוע.
    """
    data = request.json
    conditions = data.get('conditions')
    
    # קבלת רכיבי תאריך התחלה
    start_year = data.get('start_year')
    start_month = data.get('start_month')
    start_day = data.get('start_day')
    is_start_bce = data.get('is_start_bce')

    # קבלת רכיבי תאריך סיום
    end_year = data.get('end_year')
    end_month = data.get('end_month')
    end_day = data.get('end_day')
    is_end_bce = data.get('is_end_bce')

    resolution = data.get('resolution')

    # **חדש/משופר: בדיקות קלט חזקות יותר**
    if any(val is None for val in [conditions, start_year, start_month, start_day, is_start_bce,
                                   end_year, end_month, end_day, is_end_bce, resolution]):
        return jsonify({"error": "נתונים חסרים: ודא שכל שדות התנאים, התאריכים והרזולוציה סופקו."}), 400

    try:
        start_year_astro = convert_year_to_astronomical(int(start_year), bool(is_start_bce))
        end_year_astro = convert_year_to_astronomical(int(end_year), bool(is_end_bce))

        start_time = time.time()
        result = find_complex_constellation(
            conditions, 
            start_year_astro, int(start_month), int(start_day),
            end_year_astro, int(end_month), int(end_day),
            resolution
        )
        end_time = time.time()
        duration_ms = round((end_time - start_time) * 1000, 2)

        response_data = {'results': result, 'duration_ms': duration_ms}
        return jsonify(response_data)
    except ValueError as e:
        return jsonify({'error': f"שגיאה בעיבוד תאריך או תנאים: {e}"}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500


# נקודת קצה חדשה לחישוב נתונים לתרשים סינוסי/גלי (קריאה לפונקציה generate_sine_chart_data)
@app.route('/generate_sine_chart_data', methods=['POST'])
def generate_sine_chart_data_route():
    """
    מקבל שמות כוכבים ותאריכי התחלה/סיום (שנה, חודש, יום, ודגל is_bce)
    ומייצר נתונים עבור תרשים סינוסי. מחזיר את הנתונים בפורמט JSON, כולל זמן ביצוע.
    """
    data = request.json
    planet_names = data.get('planet_names')
    
    # קבלת רכיבי תאריך התחלה
    start_year = data.get('start_year')
    start_month = data.get('start_month')
    start_day = data.get('start_day')
    is_start_bce = data.get('is_start_bce')

    # קבלת רכיבי תאריך סיום
    end_year = data.get('end_year')
    end_month = data.get('end_month')
    end_day = data.get('end_day')
    is_end_bce = data.get('is_end_bce')

    aspects_to_find = data.get('aspects_to_find', [])
    resolution = data.get('resolution')
    patterns_to_find = data.get('patterns_to_find', [])

    # **חדש/משופר: בדיקות קלט חזקות יותר**
    if any(val is None for val in [planet_names, start_year, start_month, start_day, is_start_bce,
                                   end_year, end_month, end_day, is_end_bce, resolution]):
        return jsonify({"error": "נתונים חסרים: ודא שכל שדות הכוכבים, התאריכים והרזולוציה סופקו."}), 400

    try:
        start_year_astro = convert_year_to_astronomical(int(start_year), bool(is_start_bce))
        end_year_astro = convert_year_to_astronomical(int(end_year), bool(is_end_bce))

        start_time = time.time()
        result = generate_sine_chart_data(
            planet_names, 
            start_year_astro, int(start_month), int(start_day),
            end_year_astro, int(end_month), int(end_day),
            aspects_to_find, resolution, patterns_to_find
        )
        end_time = time.time()
        duration_ms = round((end_time - start_time) * 1000, 2)

        response_data = {'data_points': result, 'duration_ms': duration_ms}
        return jsonify(response_data)
    except ValueError as e:
        return jsonify({'error': f"שגיאה בעיבוד תאריך: {e}"}), 400
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# נקודת קצה לרישום תוצאות בדיקה (לצורך ניפוי באגים)
@app.route('/log_test_result', methods=['POST'])
def log_test_result():
    test_data = request.json
    timestamp = datetime.now().isoformat()
    log_entry = {"timestamp": timestamp, "test_info": test_data}
    
    # Define a specific file for logging test results
    log_file_path = 'test_results.json'

    # Check if the file exists and has content
    try:
        with open(log_file_path, 'r', encoding='utf-8') as f:
            file_content = f.read()
            if file_content.strip():
                test_results = json.loads(file_content)
            else:
                test_results = []
    except FileNotFoundError:
        test_results = []
    except json.JSONDecodeError:
        # Handle case where file is empty or contains invalid JSON
        test_results = []

    test_results.append(log_entry)

    try:
        with open(log_file_path, 'w', encoding='utf-8') as f:
            json.dump(test_results, f, indent=2, ensure_ascii=False)
        return jsonify({"status": "success"}), 200
    except Exception as e:
        app.logger.error(f"Failed to log test result to file: {e}")
        return jsonify({"status": "error", "message": str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)

