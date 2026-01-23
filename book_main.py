import os
import json
import smtplib
import ssl
import urllib.request
import urllib.error
import random
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def get_book_recommendation():
    """
    urllib을 사용하여 Google Gemini API (REST)를 직접 호출합니다.
    """
    # API 키 공백 제거
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    
    if not api_key:
        return "오류: GEMINI_API_KEY가 설정되지 않았습니다."

    # 1. 원장님의 관심사 테마 리스트 (반영 완료)
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

    # 2. 프롬프트 구성
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

    # [중요 변경] 모델명을 'gemini-2.0-flash'로 변경 (1.5 버전 종료됨)
    model_name = "gemini-2.0-flash"
    
    # 3. API 엔드포인트 설정 (v1beta + 2.0 Flash)
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"

    # 4. 요청 데이터 구성
    payload = {
        "contents": [{
            "parts": [{"text": prompt}]
        }]
    }
    data = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}

    # 5. API 호출
    try:
        req = urllib.request.Request(url, data=data, headers=headers, method="POST")
        with urllib.request.urlopen(req) as response:
            response_body = response.read().decode("utf-8")
            response_json = json.loads(response_body)
            
            try:
                text = response_json["candidates"][0]["content"]["parts"][0]["text"]
                result_text = f"Selected Theme: [{today_theme}]\n\n{text}"
                return result_text
            except (KeyError, IndexError) as e:
                return f"API 응답 파싱 실패: {e}\n응답: {response_body}"
                
    except urllib.error.HTTPError as e:
        error_content = e.read().decode("utf-8")
        return f"HTTP 에러 ({e.code}): {e.reason}\n상세: {error_content}"
    except urllib.error.URLError as e:
        return f"연결 실패: {e.reason}"
    except Exception as e:
        return f"알 수 없는 오류: {str(e)}"

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()
    
    # 받는 사람도 본인 이메일로 설정
    receiver_email = sender_email

    if not sender_email or not sender_password:
        print("이메일 설정이 누락되었습니다.")
        return

    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = receiver_email
    
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 정형외과 의사 추천 도서"

    body = f"""
    원장님, 좋은 아침입니다.
    오늘의 영감을 위한 책 추천입니다.
    
    ==================================================
    {content}
    ==================================================
    
    오늘도 평온하고 의미 있는 하루 보내시길 바랍니다.
    """
    msg.attach(MIMEText(body, "plain", "utf-8"))

    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, receiver_email, msg.as_string())
        print("이메일 발송 완료")
    except Exception as e:
        print(f"이메일 발송 실패: {e}")

if __name__ == "__main__":
    print("Gemini 2.0 모델로 책 추천 생성 중...")
    recommendation = get_book_recommendation()
    
    print(recommendation)
    send_email(recommendation)
