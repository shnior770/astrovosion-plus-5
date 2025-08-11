# השתמש בתמונת בסיס של Python עם גרסה יציבה (3.10 היא בחירה טובה)
FROM python:3.10-slim-buster

# הגדר את ספריית העבודה בתוך הקונטיינר
WORKDIR /app

# העתק את קובץ ה-requirements.txt לתוך הקונטיינר
COPY requirements.txt .

# התקן את התלויות של Python (כולל pyswisseph ו-gunicorn)
# נתקין גם את כלי הבנייה החיוניים ו-libgfortran ש-pyswisseph דורש
# (הערה: ה-apt-packages.txt שקול ל-RUN apt-get install כאן)
RUN apt-get update && apt-get install -y \
    build-essential \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir -r requirements.txt

# העתק את כל שאר קבצי הפרויקט לתוך הקונטיינר
COPY . .

# העתק את תיקיית swisseph_data/ephe למיקום הנכון אם היא לא ברמה העליונה
# אם תיקיית הנתונים כבר נקראת 'ephe' ונמצאת בשורש הפרויקט
# אז ה-COPY . . כבר יעתיק אותה, ואין צורך בשורה נפרדת
# אם קבצי ה-ephe נמצאים בתוך swisseph_data, שצריך להיות בשם ephe:
# ודא שהתיקייה המקומית נקראת ephe והקבצים נמצאים ישירות בתוכה
# או אם הקבצים הם ב-swisseph_data/ephe:
# COPY swisseph_data/ephe ephe

# וודא שפקודת ה-set_ephe_path בקוד שלך עדיין היא se.set_ephe_path('ephe')
# (מכיוון שתיקיית 'ephe' תהיה בתוך '/app/ephe')

# פקודת ההפעלה של היישום באמצעות Gunicorn
CMD ["gunicorn", "app:app", "--bind", "0.0.0.0:$PORT"]
