# השתמש בתמונת בסיס קטנה של פייתון על בסיס Debian Bullseye (נתמך)
FROM python:3.11-slim-bullseye

# הגדרת משתני סביבה לאיתור קבצי האפמריס (ephe)
ENV SWISSEPH_PATH_ABS /app/ephe

# הגדרת תיקיית העבודה בתוך הקונטיינר
WORKDIR /app

# העתק את קבצי האפמריס (ephe) לתוך הקונטיינר
# וודא שתיקיית ephe קיימת בתיקיית הפרויקט המקומית שלך
COPY ephe/ /app/ephe/

# התקנת כלי בנייה ותלויות מערכת נדרשות על ידי pyswisseph וחבילות אחרות
# --no-install-recommends מפחית את גודל התמונה
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

# העתק את קבצי הדרישות והתקן אותם
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# העתק את כל שאר קבצי היישום לתוך הקונטיינר
COPY . .

# חשוף את פורט 8080 (הפורט שבו Flask ירוץ ב-Cloud Run)
EXPOSE 8080

# הגדר את פקודת ההפעלה של היישום
CMD ["python", "app.py"]
