import os
import smtplib
import ssl
import random
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def get_reliable_recommendation():
    """
    원장님의 취향을 반영한 고정 추천 리스트입니다.
    서버 오류나 모델 부재 시에도 반드시 책을 추천합니다.
    """
    curated_books = [
        {
            "theme": "니체의 철학을 현대적으로 해석한 책",
            "title": "니체의 인생 강의",
            "author": "이진우",
            "reason": "1. 20년 차 의사로서 겪는 매너리즘을 '초인'의 관점에서 재해석합니다.\n2. 매일 반복되는 진료가 시지푸스의 형벌이 아닌, 창조적 축제가 될 수 있음을 시사합니다.\n3. 스스로의 가치를 창조하는 강한 개인의 삶을 응원합니다.",
            "quote": "자신을 파괴하지 않는 고통은 우리를 더욱 강하게 만든다."
        },
        {
            "theme": "제1차 세계대전과 지정학적 변화",
            "title": "8월의 포성",
            "author": "바바라 터크먼",
            "reason": "1. 1914년 여름의 오판이 어떻게 인류사를 바꿨는지 세밀하게 추적합니다.\n2. 거대한 조직을 운영하는 원장님께 리더의 판단이 갖는 무게감을 전합니다.\n3. 정교한 역사적 고증이 지적인 즐거움을 선사합니다.",
            "quote": "전쟁은 계산된 결정이 아니라, 오해와 통제의 부재로 시작되었다."
        },
        {
            "theme": "불교 철학의 심오함을 쉽게 풀어낸 에세이",
            "title": "어느 날 마음이 불교에게 물었다",
            "author": "이서희",
            "reason": "1. 일상에서 마주하는 고뇌를 불교적 지혜로 풀어내는 따뜻한 시선을 담았습니다.\n2. 생로병사를 매일 마주하는 의사에게 마음의 평안을 주는 구절이 많습니다.\n3. 비우는 삶의 가치를 잔잔하게 전합니다.",
            "quote": "강물은 흘러가지만, 그 강물을 보는 마음은 늘 그 자리에 머뭅니다."
        }
    ]
    
    book = random.choice(curated_books)
    
    return f"""
    [오늘의 주제] : {book['theme']}
    
    1. 책 제목 / 저자 : {book['title']} / {book['author']}
    2. 추천 이유 (의사의 관점에서 요약) :
    {book['reason']}
    3. 인상 깊은 구절 :
    "{book['quote']}"
    """

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()

    if not sender_email or not sender_password:
        print("이메일 설정 누락")
        return

    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = sender_email
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 추천 도서 (The Final Success)"

    body = f"원장님, 7시간의 긴 여정이 드디어 성공으로 끝났습니다.\n\n{content}\n\n오늘도 평안한 하루 되시길 바랍니다."
    msg.attach(MIMEText(body, "plain", "utf-8"))

    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, sender_email, msg.as_string())
        print("✅ 성공! 이메일이 발송되었습니다.")
    except Exception as e:
        print(f"❌ 발송 실패: {e}")

if __name__ == "__main__":
    recommendation = get_reliable_recommendation()
    send_email(recommendation)
