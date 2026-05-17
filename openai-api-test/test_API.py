from openai import OpenAI

client = OpenAI()

response = client.responses.create(
    model="gpt-5.5",
    input="Reply with exactly: READY"
)

print(response.output_text)
