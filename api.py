# api.py
from fastapi import FastAPI, HTTPException
from datetime import date
from typing import List, Optional

# ייבוא הפונקציות הקיימות שלנו
from astro import calculate_birth_chart
from historical_pattern_finder import find_historical_pattern, find_historical_aspect, PLANETS_HEBREW_MAP, ASPECTS_HEBREW_MAP
from sine_graph_generator import generate_sine_graph_data
from library_manager import save_search, get_all_searches, search_library, delete_search

app = FastAPI(
    title="AstroViz API",
    description="API for astrological calculations and historical data analysis.",
    version="1.0.0",
)

app = FastAPI(
    title="AstroViz API",
    description="API for astrological calculations and historical data analysis.",
    version="1.0.0",
)


@app.post("/birth_chart")
async def api_calculate_birth_chart(
    year: int, month: int, day: int, lat: float, lon: float
):
    """
    מחשב ומחזיר מפת לידה מלאה (כוכבים ובתים) עבור תאריך ומיקום נתונים.
    """
    try:
        target_date = date(year, month, day)
        chart_data = calculate_birth_chart(target_date, lat, lon)
        return chart_data
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid date: {e}")


@app.post("/historical_pattern")
async def api_find_historical_pattern(
    planet_name: str, sign_name: str, start_year: int, end_year: int
):
    """
    איתור מופעים היסטוריים של כוכב במזל מסוים בטווח תאריכים.
    """
    try:
        start_date = date(start_year, 1, 1)
        end_date = date(end_year, 12, 31)
        if planet_name not in PLANETS_HEBREW_MAP.keys():
            raise HTTPException(status_code=400, detail="Invalid planet name.")
        results = find_historical_pattern(planet_name, sign_name, start_date, end_date)
        return results
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid date: {e}")


@app.post("/historical_aspect")
async def api_find_historical_aspect(
    planet1: str, planet2: str, aspect_name: str, start_year: int, end_year: int
):
    """
    איתור מופעים היסטוריים של היבט בין שני כוכבים בטווח תאריכים.
    """
    try:
        start_date = date(start_year, 1, 1)
        end_date = date(end_year, 12, 31)
        if aspect_name not in ASPECTS_HEBREW_MAP.keys():
            raise HTTPException(status_code=400, detail="Invalid aspect name.")
        results = find_historical_aspect(planet1, planet2, aspect_name, start_date, end_date)
        return results
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid date: {e}")


@app.post("/sine_graph")
async def api_generate_sine_graph_data(
    planet_names: List[str], start_year: int, end_year: int, aspects_to_find: Optional[List[str]] = None
):
    """
    מייצר נתוני מיקום והיבטים עבור גרף סינוסי.
    """
    try:
        start_date = date(start_year, 1, 1)
        end_date = date(end_year, 12, 31)
        graph_data = generate_sine_graph_data(planet_names, start_date, end_date, aspects_to_find)
        return graph_data
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Invalid date or parameters: {e}")


@app.get("/library/all")
async def api_get_all_searches():
    """
    שולף את כל החיפושים השמורים מהספרייה.
    """
    return get_all_searches()


@app.post("/library/search")
async def api_search_library(
    query: Optional[str] = None, by_tool: Optional[str] = None, by_tag: Optional[str] = None
):
    """
    מבצע חיפוש בספרייה על בסיס מילות מפתח, כלי או תגית.
    """
    return search_library(query, by_tool, by_tag)


@app.post("/library/save")
async def api_save_search(tool_name: str, title: str, search_params: dict, results: dict, tags: Optional[List[str]] = None):
    """
    שומר חיפוש חדש בספרייה.
    """
    save_search(tool_name, title, search_params, results, tags)
    return {"message": "Search saved successfully."}


@app.delete("/library/delete/{search_id}")
async def api_delete_search(search_id: int):
    """
    מוחק חיפוש מהספרייה על פי ה-ID שלו.
    """
    delete_search(search_id)
    return {"message": f"Search with ID {search_id} deleted successfully."}
