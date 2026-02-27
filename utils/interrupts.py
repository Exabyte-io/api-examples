import asyncio
import inspect
import uuid
from dataclasses import dataclass
from functools import wraps
from typing import Any, Awaitable, Callable

from mat3ra.utils.jupyterlite.environment import ENVIRONMENT, EnvironmentsEnum

try:
    from IPython.display import HTML, display  # type: ignore
except Exception:
    HTML = None
    display = None


class UserAbortError(RuntimeError):
    pass


def display_abort_controls_in_current_cell_output(
    channel_name: str = "mat3ra_abort_channel",
    abort_button_text: str = "Abort",
) -> None:
    """
    Shows:
      [Abort]  Press ESC to abort

    Only works in notebook frontends that support HTML output.
    Safe no-op otherwise.
    """
    if HTML is None or display is None:
        return

    element_id = f"abort_controls_{uuid.uuid4().hex}"

    display(
        HTML(
            f"""
            <div style="display:flex; align-items:center; gap:12px; margin:8px 0;">
              <button
                id="{element_id}_button"
                style="
                  background:#d32f2f; color:white; border:none; padding:8px 14px;
                  border-radius:6px; cursor:pointer; font-weight:600;
                "
              >{abort_button_text}</button>

              <span style="font-family:monospace; opacity:0.85;">Press ESC to abort</span>
              <span id="{element_id}_status" style="font-family:monospace; opacity:0.85;"></span>
            </div>

            <script>
            (function() {{
              const channelName = {channel_name!r};

              // Install ESC broadcaster once per page
              if (!window.__mat3raEscapeAbortInstalled) {{
                window.__mat3raEscapeAbortInstalled = true;
                const escChannel = new BroadcastChannel(channelName);
                document.addEventListener("keydown", (event) => {{
                  if (event.key === "Escape") {{
                    escChannel.postMessage({{ type: "abort", source: "escape" }});
                  }}
                }}, true);
              }}

              // Button broadcaster (this output)
              const buttonChannel = new BroadcastChannel(channelName);
              const buttonElement = document.getElementById("{element_id}_button");
              const statusElement = document.getElementById("{element_id}_status");
              if (!buttonElement) return;

              buttonElement.addEventListener("click", () => {{
                buttonChannel.postMessage({{ type: "abort", source: "button" }});
                if (statusElement) {{
                  statusElement.textContent = "Abort sent";
                  statusElement.style.color = "#d32f2f";
                }}
              }});
            }})();
            </script>
            """
        )
    )


@dataclass
class BroadcastChannelAbortController:
    """
    WebWorker-side receiver. Works only in pyodide (emscripten).
    In regular Python: start() does nothing and is_aborted stays False.
    """

    channel_name: str = "mat3ra_abort_channel"
    is_aborted: bool = False

    def __post_init__(self) -> None:
        self._broadcast_channel = None
        self._on_message_proxy = None

    def start(self) -> None:
        if ENVIRONMENT != EnvironmentsEnum.PYODIDE:
            return
        if self._broadcast_channel is not None:
            return

        import js  # type: ignore
        from pyodide.ffi import create_proxy  # type: ignore

        self._broadcast_channel = js.BroadcastChannel.new(self.channel_name)

        def on_message(event) -> None:
            message = getattr(event, "data", None)
            if message and getattr(message, "type", None) == "abort":
                self.is_aborted = True

        self._on_message_proxy = create_proxy(on_message)
        self._broadcast_channel.onmessage = self._on_message_proxy  # type: ignore

    def stop(self) -> None:
        if self._broadcast_channel is None:
            return

        self._broadcast_channel.close()
        self._broadcast_channel = None

        if self._on_message_proxy is not None:
            self._on_message_proxy.destroy()
            self._on_message_proxy = None


async def run_interruptible_loop_async(
    loop_body: Callable[[], Awaitable[bool]],
    poll_interval_seconds: float,
    *,
    channel_name: str = "mat3ra_abort_channel",
    check_interval_seconds: float = 0.05,
    show_controls: bool = True,
) -> None:
    """
    Minimal wrapper.

    loop_body():
      - do one "poll" iteration
      - return True to keep looping, False to stop normally

    Between iterations we sleep in small slices so:
      - pyodide: ESC/button can be received and stop the loop
      - regular Python: yields control (Ctrl+C/Stop works where supported)
    """
    broadcast_channel_abort_controller = BroadcastChannelAbortController(channel_name=channel_name)
    broadcast_channel_abort_controller.start()

    if show_controls and ENVIRONMENT != EnvironmentsEnum.PYODIDE:
        display_abort_controls_in_current_cell_output(channel_name=channel_name, abort_button_text="Abort")

    try:
        while True:
            should_continue = await loop_body()
            if not should_continue:
                return

            remaining_seconds = float(poll_interval_seconds)
            while remaining_seconds > 0:
                if broadcast_channel_abort_controller.is_aborted:
                    raise UserAbortError("Aborted by user.")
                await asyncio.sleep(min(check_interval_seconds, remaining_seconds))
                remaining_seconds -= check_interval_seconds

    finally:
        broadcast_channel_abort_controller.stop()


def interruptible_polling_loop(
    poll_interval_kwarg_name: str = "poll_interval",
    *,
    default_poll_interval_seconds: float = 10.0,
    channel_name: str = "mat3ra_abort_channel",
    check_interval_seconds: float = 0.05,
    show_controls: bool = True,
):
    """
    Decorator for single-iteration polling functions.

    The decorated function must:
      - return True to continue
      - return False to stop

    The decorated function becomes an `async def` that runs the polling loop until completion.

    The polling interval can be passed at call time using `poll_interval_kwarg_name` (defaults to
    `"poll_interval"`). If not provided, `default_poll_interval_seconds` is used.
    """

    def decorator(poll_step_function: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(poll_step_function)
        async def wrapped(*args: Any, **kwargs: Any) -> None:
            poll_interval_seconds = float(kwargs.pop(poll_interval_kwarg_name, default_poll_interval_seconds))

            async def loop_body() -> bool:
                result = poll_step_function(*args, **kwargs)
                should_continue = await result if inspect.isawaitable(result) else result
                return bool(should_continue)

            await run_interruptible_loop_async(
                loop_body,
                poll_interval_seconds,
                channel_name=channel_name,
                check_interval_seconds=check_interval_seconds,
                show_controls=show_controls,
            )

        return wrapped

    return decorator
