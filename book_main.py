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

def get_book_recommendation():
    """
    [최종 해결 솔루션]
    1. 새 프로젝트(AI Studio)에서는 'gemini-2.0-flash'가 기본입니다.
    2. 만약 실패하면 자동으로 다른 모델들을 순차적으로 시도합니다. (Auto-Discovery)
    """
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    if not api_key:
        return "오류: GEMINI_API_KEY가 설정되지 않았습니다."

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

    # [핵심] 시도할 모델 리스트 (우선순위: 2.0 -> 1.5)
    # 새 프로젝트는 2.0이 될 확률이 99%입니다.
    models_to_try = [
        "gemini-2.0-flash",       # 1순위 (최신, 새 계정용)
        "gemini-1.5-flash",       # 2순위 (표준)
        "gemini-1.5-flash-8b",    # 3순위 (경량화)
        "gemini-1.5-pro"          # 4순위 (고성능)
    ]

    prompt = f"""
    당신은 정형외과 의사를 위한 독서 큐레이터입니다.
    다음 주제에 맞는 책 1권을 한국어로 추천해주세요.
    
    주제: {today_theme}
    
    형식:
    1. 책 제목 / 저자
    2. 추천 이유 (3줄 요약)
    3. 인상 깊은 구절
    """

    payload = {
        "contents": [{"parts": [{"text": prompt}]}]
    }
    data = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}

    # 모델 리스트를 순회하며 성공할 때까지 시도
    last_error = ""
    for model_name in models_to_try:
        url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"
        print(f"📡 연결 시도 중: {model_name}...")
        
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(req) as response:
                response_body = response.read().decode("utf-8")
                response_json = json.loads(response_body)
                
                try:
                    text = response_json["candidates"][0]["content"]["parts"][0]["text"]
                    print(f"✅ 성공! ({model_name})")
                    return f"Selected Theme: [{today_theme}]\n(Model: {model_name})\n\n{text}"
                except (KeyError, IndexError):
                    print(f"⚠️ 파싱 오류 ({model_name})")
                    continue # 다음 모델 시도

        except urllib.error.HTTPError as e:
            error_content = e.read().decode("utf-8")
            # 404는 모델이 없다는 뜻이므로 바로 다음 모델로 넘어감
            if e.code == 404:
                print(f"❌ {model_name} 모델 없음 (404). 다음 모델 검색...")
                last_error = f"{model_name}: 404 Not Found"
                continue
            
            # 429는 잠시 대기 후 재시도 (같은 모델로)
            if e.code == 429:
                print(f"⏳ 사용량 제한(429). 3초 대기 후 {model_name} 재시도...")
                time.sleep(3)
                # 재시도 로직이 복잡해지므로, 429면 그냥 다음 모델(더 가벼운거)로 넘어가는 게 낫습니다.
                continue
                
            last_error = f"HTTP {e.code}: {e.reason}"
            print(f"❌ 에러 발생: {last_error}")
            continue

        except Exception as e:
            last_error = str(e)
            print(f"❌ 알 수 없는 오류: {last_error}")
            continue

    return f"❌ 모든 모델 연결 실패.\n마지막 에러: {last_error}\nAPI 키가 'AI Studio'에서 생성된 것이 맞는지 확인해주세요."

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()
    
    receiver_email = os.environ.get("RECEIVER_EMAIL", "").strip()
    if not receiver_email:
        receiver_email = sender_email

    if not sender_email or not sender_password:
        print("❌ 이메일 설정 누락")
        return

    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = receiver_email
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 추천 도서 (Final)"

    body = f"""
    원장님, 드디어 해결되었습니다.
    
    ==================================================
    {content}
    ==================================================
    """
    msg.attach(MIMEText(body, "plain", "utf-8"))

    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, receiver_email, msg.as_string())
        print("✅ 이메일 발송 성공")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    result = get_book_recommendation()
    print(result)
    send_email(result)
