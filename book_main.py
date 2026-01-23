import smtplib
import os
import random
from email.mime.text import MIMEText
from email.header import Header
from datetime import datetime
import google.generativeai as genai

# ==========================================
# [설정] 환경 변수
# ==========================================
MY_EMAIL = os.getenv("MY_EMAIL")
MY_APP_PASSWORD = os.getenv("MY_APP_PASSWORD")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")

# 받는 사람 (없으면 나에게)
env_receiver = os.getenv("RECEIVER_EMAIL")
RECEIVER_EMAIL = env_receiver if env_receiver else MY_EMAIL

# ==========================================
# [설정] Gemini 연결 (선생님 픽: 2.0 모델)
# ==========================================
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY.strip())

def get_book_recommendation():
    if not GEMINI_API_KEY:
        return "<h3>설정 오류</h3><p>API 키가 없습니다.</p>"

    themes = [
        "니체의 철학을 현대적으로 해석한 책",
        "제1차 세계대전과 지정학적 변화를 다룬 역사서",
        "불교 철학의 심오함을 쉽게 풀어낸 에세이",
        "빅 히스토리(Big History) 관점의 인문학",
        "실존주의 문학의 숨겨진 명작",
        "인공지능과 인간의 미래를 다룬 통찰력 있는 책",
        "일본의 미학(Wabi-sabi)이나 장인 정신을 다룬 책",
        "의사가 쓴 죽음과 삶에 대한 에세이"
    ]
    
    today_theme = random.choice(themes)

    prompt = f"""
    당신은 지적인 50대 정형외과 의사를 위한 '독서 큐레이터'입니다.
    
    [오늘의 주제] : {today_theme}
    
    위 주제와 관련하여, 깊이 있고 통찰력을 주는 책 1권을 추천해주세요.
    (너무 뻔한 베스트셀러보다는, 숨겨진 명작이나 깊이 있는 책을 선호합니다.)

    [출력 형식]
    1. 책 제목 / 저자
    2. 추천 이유 (의사의 관점에서 흥미로울 포인트 3줄 요약)
    3. 인상 깊은 구절 (한 문장)
    """

    try:
        # [핵심] 아까 연결은 성공했던 그 모델 (2.0)
        model = genai.GenerativeModel('gemini-2.0-flash')
        response = model.generate_content(prompt)
        return f"<h3>🎨 오늘의 테마: {today_theme}</h3><hr><br>" + response.text.replace('\n', '<br>')
    except Exception as e:
