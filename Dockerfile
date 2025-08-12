# השתמש בתמונת בסיס של Python עם גרסה יציבה (3.10 היא בחירה טובה)
FROM python:3.10-slim-bullseye

# הגדר את ספריית העבודה בתוך הקונטיינר
WORKDIR /app

# העתק את קובץ ה-requirements.txt לתוך הקונטיינר
COPY requirements.txt .

# התקן את התלויות של Python (כולל pyswisseph ו-gunicorn)
# נתקין גם את כלי הבנייה החיוניים ו-libgfortran ש-pyswisseph דורש
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir -r requirements.txt \
    && pip show pyswisseph  # הוסף זאת כדי לוודא התקנה

# **בדיקה מפורשת של הייבוא במהלך הבנייה**
# אם זה נכשל, ה-Build כולו ייכשל כאן וזה ייתן לנו אבחנה מדויקת יותר


# העתק את כל שאר קבצי הפרויקט לתוך הקונטיינר
COPY . .

# פקודת ההפעלה של היישום באמצעות Gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "app:app"]
