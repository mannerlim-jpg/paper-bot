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
    모델: HuggingFaceH4/zephyr-7b-beta (무료 서버에서 가장 안정적이고 성능 좋은 모델)
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

    # 2. Hugging Face API 주소 및 모델 설정 (Zephyr 모델 사용)
    model_id = "HuggingFaceH4/zephyr-7b-beta"
    url = f"https://api-inference.huggingface.co/models/{model_id}"
    
    # 3. 프롬프트 (Zephyr 모델은 <|system|> 태그를 좋아합니다)
    prompt = f"""<|system|>
    당신은 20년 차 정형외과 의사를 위한 지적인 '독서 큐레이터'입니다.
    반드시 '한국어'로 답변해야 합니다.</s>
    <|user|>
    아래 주제에 맞춰 깊이 있는 책 1권을 추천해주세요.

    주제: {today_theme}

    [필수 포함 내용]
    1. 책 제목과 저자
    2. 추천 이유 (의사의 관점에서 3줄 요약)
    3. 인상 깊은 구절 (한 문장)</s>
    <|assistant|>
    """

    payload = {
        "inputs": prompt,
        "parameters": {
            "max_new_tokens": 1000,
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
                
                # [수정 완료] 아까 잘렸던 부분이 여기입니다!
                response_json = json.loads(response_body)
                
                if isinstance(response_json, list) and "generated_text" in response_json[0]:
                    text = response_json[0]["generated_text"]
                    return f"Selected Theme: [{today_theme}]\n\n{text}"
                else:
                    return f"응답 형식 이상함: {response_body}"

        except urllib.error.HTTPError as e:
            error_msg = e.read().decode("utf-8")
            # 503 (로딩중) 처리 - 무료 모델은 깨우는데 시간이 걸릴 수 있음
            if e.code == 503:
                print(f"⚠️ 모델 로딩 중 (503)... 20초 후 재시도 ({attempt+1}/{max_retries})")
                time.sleep(20)
                continue
            return f"HTTP 에러 ({e.code}): {error_msg}"
        except Exception as e:
            return f"알 수 없는 오류: {str(e)}"

    return "❌ 서버 연결 실패 (잠시 후 다시 시도해주세요)"

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()
    receiver_email = sender_email

    if not sender_email or not sender_password:
        print("❌ [오류] 이메일 설정이 누락되었습니다.")
        return

    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = receiver_email
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 추천 도서 (Zephyr AI)"

    body = f"""
    원장님, 좋은 아침입니다.
    오늘의 영감을 위한 책 추천입니다. (Powered by Zephyr)
    
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
        print("✅ 이메일 발송 성공!")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    print("책 추천 생성 시작...")
    recommendation = get_book_recommendation_huggingface()
    print("결과 내용:\n", recommendation)
    send_email(recommendation)
