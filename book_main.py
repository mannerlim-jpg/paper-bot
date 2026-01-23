import os
import json
import smtplib
import ssl
import urllib.request
import urllib.error
import random
import time
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def get_book_recommendation_huggingface():
    """
    Hugging Face의 무료 Inference API를 사용합니다.
    URL이 최신 주소(router.huggingface.co)로 업데이트되었습니다.
    """
    # GitHub Secrets에서 가져온 토큰
    api_token = os.environ.get("HF_TOKEN", "").strip()
    
    if not api_token:
        print("❌ [오류] 시스템 환경변수에서 HF_TOKEN을 찾을 수 없습니다.")
        return "오류: HF_TOKEN이 설정되지 않았습니다. (YAML 파일을 확인해주세요)"

    # 1. 테마 선택
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

    # 2. Hugging Face API 주소 (최신 주소로 변경됨)
    model_id = "mistralai/Mistral-7B-Instruct-v0.3"
    # [수정] api-inference -> router 로 변경
    url = f"https://router.huggingface.co/models/{model_id}"
    
    # 3. 프롬프트 (한국어 답변 요청)
    prompt = f"""<s>[INST] 당신은 20년 차 정형외과 의사를 위한 지적인 '독서 큐레이터'입니다.
    아래 주제에 맞춰 깊이 있는 책 1권을 추천해주세요. 
    반드시 '한국어'로 답변해야 합니다.

    주제: {today_theme}

    [필수 포함 내용]
    1. 책 제목과 저자
    2. 추천 이유 (의사의 관점에서 3줄 요약)
    3. 인상 깊은 구절 (한 문장) [/INST]
    """

    payload = {
        "inputs": prompt,
        "parameters": {
            "max_new_tokens": 1000, # 답변 길이 조금 더 넉넉하게
            "temperature": 0.7,
            "return_full_text": False
        }
    }
    data = json.dumps(payload).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_token}"
    }

    # 4. 재시도 로직
    max_retries = 5
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(req) as response:
                response_body = response.read().decode("utf-8")
                response_json =
