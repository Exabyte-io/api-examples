import asyncio

import pytest
from mat3ra.notebooks_utils.pyodide.runtime import (
    UserAbortError,
    interruptible_polling_loop,
    run_interruptible_loop_async,
)

POLL_INTERVAL_SECONDS = 0.01
CHECK_INTERVAL_SECONDS = 0.005


@pytest.mark.asyncio
async def test_run_interruptible_loop_async_stops_when_body_returns_false():
    call_count = 0

    async def loop_body():
        nonlocal call_count
        call_count += 1
        return call_count < 3

    await run_interruptible_loop_async(
        loop_body,
        POLL_INTERVAL_SECONDS,
        show_controls=False,
        check_interval_seconds=CHECK_INTERVAL_SECONDS,
    )
    assert call_count == 3


@pytest.mark.asyncio
async def test_interruptible_polling_loop_decorator_returns_coroutine_and_runs_until_false():
    call_count = 0

    @interruptible_polling_loop(show_controls=False)
    def poll_step():
        nonlocal call_count
        call_count += 1
        return call_count < 2

    assert asyncio.iscoroutinefunction(poll_step)
    await poll_step(poll_interval=POLL_INTERVAL_SECONDS)
    assert call_count == 2


def test_user_abort_error_is_runtime_error():
    error = UserAbortError("test message")
    assert isinstance(error, RuntimeError)
    assert str(error) == "test message"
